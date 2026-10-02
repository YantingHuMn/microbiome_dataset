def find_candidate_files(paper_id, output_dir, data_sources=None, timeout=30):
    import json
    import re
    import shutil
    import urllib.parse
    import urllib.request
    import xml.etree.ElementTree as ET
    import zipfile
    from email.message import Message
    from html.parser import HTMLParser
    from pathlib import Path

    import pandas as pd

    def request_url(url):
        request = urllib.request.Request(
            url, headers={"User-Agent": "paper-data-source-finder/1.0"}
        )
        return urllib.request.urlopen(request, timeout=timeout)

    # 1. Find first author's last name and publication year
    pid = re.sub(
        r"^https?://(?:dx\.)?doi\.org/", "", str(paper_id).strip(), flags=re.I
    )

    if re.fullmatch(r"PMC\d+", pid, re.I):
        query = f"EXT_ID:{pid.upper()} AND SRC:PMC"
        id_field = "pmcid"
    elif pid.isdigit():
        query = f"EXT_ID:{pid} AND SRC:MED"
        id_field = "id"
    elif re.fullmatch(r"10\.\d{4,9}/\S+", pid):
        query = f'DOI:"{pid}"'
        id_field = "doi"
    else:
        raise ValueError("paper_id must be DOI, numeric PMID or PMCID")

    search_url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        + urllib.parse.urlencode({
            "query": query,
            "format": "json",
            "resultType": "core",
        })
    )

    with request_url(search_url) as response:
        records = json.loads(response.read()).get(
            "resultList", {}
        ).get("result", [])

    record = next(
        (
            item for item in records
            if str(item.get(id_field, "")).lower() == pid.lower()
        ),
        None,
    )
    if record is None:
        raise RuntimeError(f"Paper metadata not found: {paper_id}")

    authors = record.get("authorList", {}).get("author", [])
    last_name = authors[0].get("lastName", "") if authors else ""
    year = str(record.get("pubYear", ""))

    if not last_name or not year:
        raise RuntimeError("First author's last name or publication year is missing")

    def safe_name(name):
        name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", name)
        return name.strip(" .") or "unnamed"

    # 2. Reuse the directory registered for this paper
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    index_file = output_dir / "paper_index.json"

    if index_file.exists():
        with index_file.open("r", encoding="utf-8") as handle:
            paper_index = json.load(handle)
    else:
        paper_index = {}

    # DOI 优先，避免同一篇分别输入 DOI、PMID 时重复建目录
    paper_key = (
        record.get("doi")
        or record.get("pmcid")
        or pid
    ).strip().lower()

    def letter_suffix(number):
        suffix = ""
        while number:
            number, remainder = divmod(number - 1, 26)
            suffix = chr(97 + remainder) + suffix
        return suffix

    if paper_key in paper_index:
        folder_name = paper_index[paper_key]["folder_name"]
        paper_dir = output_dir / folder_name
        paper_dir.mkdir(exist_ok=True)
    else:
        base_name = safe_name(f"{last_name}{year}")
        used_names = {
            item["folder_name"] for item in paper_index.values()
        }
        index = 0

        while True:
            folder_name = (
                base_name if index == 0
                else f"{base_name}_{letter_suffix(index)}"
            )
            paper_dir = output_dir / folder_name

            if folder_name in used_names or paper_dir.exists():
                index += 1
                continue

            try:
                paper_dir.mkdir()
                break
            except FileExistsError:
                index += 1

        paper_index[paper_key] = {
            "folder_name": folder_name,
            "title": record.get("title", ""),
        }

        temporary_file = index_file.with_suffix(".json.tmp")
        with temporary_file.open("w", encoding="utf-8") as handle:
            json.dump(paper_index, handle, ensure_ascii=False, indent=2)
        temporary_file.replace(index_file)

    original_dir = paper_dir / "Original Data"
    original_dir.mkdir(exist_ok=True)

    if data_sources is None:
        data_sources = find_data_sources(paper_id, timeout=timeout)

    files = []
    issues = []
    source_by_path = {}

    file_extensions = (
        ".xlsx", ".xls", ".csv", ".tsv", ".txt",
        ".zip", ".rds", ".rdata", ".biom", ".h5ad",
        ".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa",
    )

    class FileLinkParser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.links = []

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                href = dict(attrs).get("href")
                if href:
                    self.links.append(href)

    def unique_path(directory, filename):
        path = directory / filename
        index = 1
        while path.exists():
            path = directory / f"{Path(filename).stem}_{index}{Path(filename).suffix}"
            index += 1
        return path

    # 3. Download attachments; inspect HTML for direct file links
    visited = set()

    for source in data_sources:
        source_link = source.get("link") or source.get("value")

        if source.get("type") == "raw_reads":
            files.append({
                **source,
                "local_path": None,
                "candidate_type": "raw_reads",
                "reason": "Raw reads entry recorded; sequencing files not downloaded",
            })
            continue

        if not source_link:
            issues.append({"source": source, "reason": "Missing link"})
            continue

        queue = [(source_link, 0)]

        while queue:
            url, depth = queue.pop(0)
            if url in visited:
                continue
            visited.add(url)

            try:
                with request_url(url) as response:
                    actual_url = response.geturl()
                    content_type = response.headers.get("Content-Type", "").lower()
                    head = response.read(512)

                    is_html = (
                        "text/html" in content_type
                        or head.lstrip().lower().startswith(
                            (b"<!doctype html", b"<html")
                        )
                    )

                    if is_html:
                        # when PMC supplementary return HTML, try PLOS download entrance for the same file
                        filename = urllib.parse.unquote(
                            urllib.parse.urlsplit(url).path
                        ).rsplit("/", 1)[-1]

                        article_doi = record.get("doi", "")
                        if (
                            urllib.parse.urlsplit(url).hostname
                            == "pmc.ncbi.nlm.nih.gov"
                            and article_doi.startswith("10.1371/journal.")
                        ):
                            article_name = article_doi.split(
                                "10.1371/journal.", 1
                            )[1]

                            if re.fullmatch(
                                re.escape(article_name)
                                + r"\.s\d+\.[A-Za-z0-9]+",
                                filename,
                            ):
                                supplement_id = (
                                    "10.1371/journal."
                                    + filename.rsplit(".", 1)[0]
                                )
                                fallback_url = (
                                    "https://journals.plos.org/"
                                    + (
                                        "plosone"
                                        if article_name.startswith("pone.")
                                        else article_name.split(".", 1)[0]
                                    )
                                    + "/article/file?"
                                    + urllib.parse.urlencode({
                                        "type": "supplementary",
                                        "id": supplement_id,
                                    })
                                )
                                queue.append((fallback_url, depth + 1))
                                continue

                        html = (head + response.read()).decode(
                            "utf-8", errors="replace"
                        )
                        parser = FileLinkParser()
                        parser.feed(html)

                        links = []
                        for href in parser.links:
                            link = urllib.parse.urljoin(actual_url, href)
                            parsed = urllib.parse.urlsplit(link)
                            if parsed.scheme not in {"http", "https"}:
                                continue

                            target = urllib.parse.unquote(
                                parsed.path + "?" + parsed.query
                            ).lower()

                            if (
                                any(
                                    re.search(re.escape(ext) + r"(?:$|[?&#])", target)
                                    for ext in file_extensions
                                )
                                or re.search(
                                    r"(?:/download(?:/|$)|"
                                    r"type=supplementary|/article/file\?)",
                                    target,
                                )
                            ):
                                links.append(link)

                        links = list(dict.fromkeys(links))

                        if links and depth < 2:
                            queue.extend((link, depth + 1) for link in links)
                        else:
                            issues.append({
                                "link": url,
                                "reason": (
                                    "No direct file link found, page blocked, "
                                    "or platform-specific listing rule needed"
                                ),
                            })
                        continue

                    message = Message()
                    message["Content-Disposition"] = response.headers.get(
                        "Content-Disposition", ""
                    )
                    filename = message.get_filename()

                    if not filename:
                        filename = urllib.parse.unquote(
                            urllib.parse.urlsplit(actual_url).path
                        ).rsplit("/", 1)[-1]

                    filename = safe_name(filename or "download")
                    if not Path(filename).suffix:
                        if head.startswith(b"PK\x03\x04"):
                            filename += ".zip"
                        elif "csv" in content_type:
                            filename += ".csv"

                    path = unique_path(original_dir, filename)
                    try:
                        with path.open("wb") as handle:
                            handle.write(head)
                            shutil.copyfileobj(response, handle)
                    except Exception:
                        path.unlink(missing_ok=True)
                        raise

                    source_by_path[path] = {
                        "link": url,
                        "source_link": source_link,
                        "source_section": source.get("source_section"),
                        "platform": source.get("platform"),
                    }

            except Exception as exc:
                issues.append({"link": url, "reason": str(exc)})

    # 4. Unzip archives, including nested ZIPs
    pending = list(source_by_path)
    unpacked = set()

    while pending:
        path = pending.pop(0)
        if path in unpacked or not zipfile.is_zipfile(path):
            continue

        # XLSX and some other formats are ZIP containers, not standalone archives
        if path.suffix.lower() in {".xlsx", ".docx", ".pptx"}:
            continue

        unpacked.add(path)
        extract_dir = unique_path(path.parent, f"{path.stem}_extracted")
        extract_dir.mkdir()

        try:
            with zipfile.ZipFile(path) as archive:
                for member in archive.infolist():
                    relative = Path(member.filename.replace("\\", "/"))

                    if relative.is_absolute() or ".." in relative.parts:
                        raise ValueError(f"Unsafe ZIP path: {member.filename}")

                    target = extract_dir / relative

                    if member.is_dir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue

                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(member) as source_handle:
                        with target.open("wb") as target_handle:
                            shutil.copyfileobj(source_handle, target_handle)

                    source_by_path[target] = {
                        **source_by_path[path],
                        "archive_path": str(path),
                    }
                    pending.append(target)

        except Exception as exc:
            issues.append({"local_path": str(path), "reason": str(exc)})

    # 5. Read tables and classify conservatively
    def inspect_table(frame, name):
        frame = frame.dropna(how="all").dropna(axis=1, how="all")

        if frame.empty:
            return "needs_review", "Empty table"

        preview = " ".join(
            [name]
            + [str(column) for column in frame.columns]
            + [
                str(value)
                for value in frame.head(8).to_numpy().ravel()
                if pd.notna(value)
            ]
        ).lower()

        if re.search(
            r"number of gene hits|kegg level|brite hierarchies", preview
        ):
            return "functional_summary", "KEGG/function annotation summary"

        if re.search(
            r"\bp[- _]?value\b|\badjusted p\b|\bpadj\b|"
            r"\bstandard deviation\b|\bstandard error\b",
            preview,
        ):
            return "statistical_summary", "Statistical comparison fields detected"

        if re.search(
            r"clean reads|filtered reads|q20|q30|sequencing depth", preview
        ):
            return "sequencing_summary", "Sequencing/QC summary fields detected"

        numeric_columns = []
        for column in frame.columns:
            values = frame[column].dropna().astype(str).str.strip()
            if values.empty:
                continue
            numbers = pd.to_numeric(values, errors="coerce")
            if numbers.notna().mean() >= 0.9:
                numeric_columns.append(column)

        if len(numeric_columns) < 2:
            return "needs_review", "Fewer than two mostly numeric columns"

        numeric = frame[numeric_columns].apply(
            pd.to_numeric, errors="coerce"
        )
        if (numeric < 0).any().any():
            return "needs_review", "Negative values; abundance not established"

        if re.search(
            r"abundance|otu|asv|rpkm|tpm|count.matrix|taxonomic.profile",
            preview,
        ):
            return (
                "possible_abundance",
                "Abundance-related text and multiple numeric columns; "
                "sample columns and feature identities still need confirmation",
            )

        return (
            "needs_review",
            "Numeric table found, but abundance meaning is not established",
        )

    for path, provenance in source_by_path.items():
        item = {
            "file_name": path.name,
            "local_path": str(path),
            **provenance,
            "candidate_type": "needs_review",
            "tables": [],
        }

        suffix = path.suffix.lower()
        name = path.name.lower()

        try:
            if path in unpacked:
                item["candidate_type"] = "archive"
                item["reason"] = "Archive unpacked; extracted files checked separately"

            elif name.endswith((
                ".fastq", ".fq", ".fastq.gz", ".fq.gz",
            )):
                item["candidate_type"] = "raw_reads"
                item["reason"] = "Sequencing read file"

            elif suffix in {".xlsx", ".xls"}:
                with pd.ExcelFile(path) as workbook:
                    for sheet in workbook.sheet_names:
                        frame = pd.read_excel(workbook, sheet_name=sheet)
                        kind, reason = inspect_table(frame, f"{path.name} {sheet}")
                        item["tables"].append({
                            "sheet": sheet,
                            "rows": len(frame),
                            "columns": len(frame.columns),
                            "candidate_type": kind,
                            "reason": reason,
                        })

            elif suffix in {".csv", ".tsv", ".txt"}:
                frame = pd.read_csv(
                    path,
                    sep="\t" if suffix == ".tsv" else None,
                    engine="python",
                )
                kind, reason = inspect_table(frame, path.name)
                item["tables"].append({
                    "sheet": None,
                    "rows": len(frame),
                    "columns": len(frame.columns),
                    "candidate_type": kind,
                    "reason": reason,
                })

            elif suffix in {".rds", ".rdata", ".biom", ".h5ad"}:
                item["reason"] = "Requires a format-specific reader"

            else:
                item["reason"] = "File format not handled by this version"

            if item["tables"]:
                kinds = {table["candidate_type"] for table in item["tables"]}
                item["candidate_type"] = (
                    "possible_abundance"
                    if "possible_abundance" in kinds
                    else next(iter(kinds)) if len(kinds) == 1
                    else "needs_review"
                )

        except Exception as exc:
            item["candidate_type"] = "needs_review"
            item["reason"] = f"Unable to read file: {exc}"

        files.append(item)

    return {
        "paper_id": paper_id,
        "paper_directory": str(paper_dir),
        "original_data_directory": str(original_dir),
        "files": files,
        "issues": issues,
    }