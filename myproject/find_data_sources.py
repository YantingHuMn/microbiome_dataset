import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser


def find_data_sources(paper_id, timeout=30):
    """Input DOI, PMID, PMCID, return Data Availability's data sources (links and raw reads accessions)."""

    def read_url(url):
        import time
        import urllib.error

        request = urllib.request.Request(
            url, headers={"User-Agent": "paper-data-source-finder/1.0"}
        )

        for attempt in range(3):
            try:
                with urllib.request.urlopen(
                    request, timeout=timeout
                ) as response:
                    return response.read()

            except Exception as exc:
                retryable = (
                    isinstance(exc, urllib.error.HTTPError)
                    and exc.code in {429, 500, 502, 503, 504}
                ) or (
                    isinstance(exc, (urllib.error.URLError, TimeoutError))
                    and not isinstance(exc, urllib.error.HTTPError)
                )

                if retryable and attempt < 2:
                    print(
                        f"请求失败，准备重试 {attempt + 1}/2：{url}；{exc}",
                        flush=True,
                    )
                    time.sleep(2 * (attempt + 1))
                    continue

                raise RuntimeError(
                    f"Failed to read: {url}；{exc}"
                ) from exc
                
    # 1. Locate PMC all full text
    pid = str(paper_id).strip()
    pid = re.sub(
        r"^https?://(?:dx\.)?doi\.org/", "", pid, flags=re.I
    )

    if re.fullmatch(r"PMC\d+", pid, re.I):
        pmcid = pid.upper()
    else:
        if pid.isdigit():
            query = f"EXT_ID:{pid} AND SRC:MED"
        elif re.fullmatch(r"10\.\d{4,9}/\S+", pid):
            query = f'DOI:"{pid}"'
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
        records = json.loads(read_url(search_url)).get(
            "resultList", {}
        ).get("result", [])

        if pid.isdigit():
            matching = [
                r for r in records
                if str(r.get("id", "")) == pid and r.get("source") == "MED"
            ]
        else:
            matching = [
                r for r in records
                if str(r.get("doi", "")).lower() == pid.lower()
            ]

        pmcid = next(
            (r["pmcid"] for r in matching if r.get("pmcid")), None
        )
        if not pmcid:
            raise RuntimeError(
                "Failed to find readable PMC full text; current version only supports PMC full text, "
                "cannot determine the absence of data based on this."
            )

    fulltext_url = (
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
        f"{pmcid}/fullTextXML"
    )
    root = ET.fromstring(read_url(fulltext_url))

    # 2. Only locate Data Availability part
    sections = []
    for node in root.iter():
        tag = node.tag.rsplit("}", 1)[-1]
        if tag not in {"sec", "notes", "fn"}:
            continue

        title = " ".join(
            "".join(child.itertext())
            for child in node
            if child.tag.rsplit("}", 1)[-1] in {"title", "label"}
        )
        section_type = " ".join(
            node.attrib.get(key, "")
            for key in ("sec-type", "notes-type", "fn-type")
        )

        if (
            re.search(
                r"data[\s_-]*(?:availability|accessibility)"
                r"|availability\s+of\s+data",
                title + " " + section_type,
                re.I,
            )
            or re.fullmatch(r"(?:associated[\s_-]+)?data", title.strip(), re.I)
            or re.fullmatch(r"(?:associated[\s_-]+)?data", section_type.strip(), re.I)
        ):
            sections.append(node)

    if not sections:
        raise RuntimeError(
            "Full text has been read, but no Data Availability section was identified;"
            "additional section identification rules may be needed."
        )

    def section_name(section):
        title = " ".join(
            "".join(child.itertext())
            for child in section
            if child.tag.rsplit("}", 1)[-1] == "title"
        ).strip()
        if title:
            return title

        label = " ".join(
            "".join(child.itertext())
            for child in section
            if child.tag.rsplit("}", 1)[-1] == "label"
        ).strip()
        if label:
            return label

        return next(
            (
                section.attrib[key]
                for key in ("sec-type", "notes-type", "fn-type", "id")
                if section.attrib.get(key)
            ),
            "unknown",
        )

    # 3. Extract links and raw reads accessions from the Data Availability section
    def platform_for(url):
        host = (urllib.parse.urlparse(url).hostname or "").lower()
        doi_prefix = r"^https?://(?:dx\.)?doi\.org/"

        if (
            "figshare" in host
            or re.search(doi_prefix + r"10\.6084/m9\.figshare\.", url, re.I)
        ):
            return "figshare"
        if (
            host == "zenodo.org" or host.endswith(".zenodo.org")
            or re.search(doi_prefix + r"10\.5281/zenodo\.", url, re.I)
        ):
            return "zenodo"
        if (
            host == "datadryad.org" or host.endswith(".datadryad.org")
            or re.search(doi_prefix + r"10\.5061/dryad\.", url, re.I)
        ):
            return "Dryad"
        if host == "osf.io" or host.endswith(".osf.io"):
            return "OSF"
        if host == "github.com" or host.endswith(".github.com"):
            return "GitHub"
        if host == "ncbi.nlm.nih.gov" or host.endswith(".ncbi.nlm.nih.gov"):
            return "NCBI"
        if host == "ebi.ac.uk" or host.endswith(".ebi.ac.uk"):
            return "EMBL-EBI"
        if host == "ddbj.nig.ac.jp" or host.endswith(".ddbj.nig.ac.jp"):
            return "DDBJ"
        if host == "ngdc.cncb.ac.cn":
            return "NGDC GSA"
        return "unknown"

    result = []
    seen = set()

    def raw_reads_link(accession, platform):
        if platform == "NGDC GSA":
            return f"https://ngdc.cncb.ac.cn/gsa/browse/{accession}"

        if platform == "NCBI SRA":
            if accession.startswith("PRJNA"):
                return f"https://www.ncbi.nlm.nih.gov/bioproject/{accession}/"
            return (
                "https://www.ncbi.nlm.nih.gov/sra?"
                + urllib.parse.urlencode({"term": accession})
            )

        if platform == "ENA":
            return f"https://www.ebi.ac.uk/ena/browser/view/{accession}"

        if platform == "DDBJ DRA":
            if accession.startswith("PRJDB"):
                resource = "bioproject"
            else:
                resource = {
                    "DRP": "sra-study",
                    "DRR": "sra-run",
                    "DRX": "sra-experiment",
                    "DRS": "sra-sample",
                }[accession[:3]]

            return f"https://ddbj.nig.ac.jp/resource/{resource}/{accession}"

        raise ValueError(f"No raw reads link rule for platform: {platform}")


    def add(kind, platform, value, source_section):
        accession = None
        link = value

        if kind == "raw_reads":
            match = re.search(
                accession_pattern, urllib.parse.unquote(value), re.I
            )
            if not match:
                raise ValueError(f"Raw reads accession not identified: {value}")

            accession = match.group().upper()

            if not re.match(r"https?://", value, re.I):
                link = raw_reads_link(accession, platform)

        key = (kind, platform, accession if kind == "raw_reads" else link)

        if key not in seen:
            seen.add(key)

            item = {
                "type": kind,
                "source_section": source_section,
                "platform": platform,
            }
            if kind == "raw_reads":
                item["accession"] = accession

            item["link"] = link
            result.append(item)

    # Follow internal references to supporting/supplementary sections
    nodes_by_id = {
        node.attrib["id"]: node
        for node in root.iter()
        if node.attrib.get("id")
    }
    selected = set(sections)

    for section in list(sections):
        for node in section.iter():
            if node.tag.rsplit("}", 1)[-1] != "xref":
                continue

            for rid in node.attrib.get("rid", "").split():
                target = nodes_by_id.get(rid)
                if target is None:
                    continue

                target_title = " ".join(
                    "".join(child.itertext())
                    for child in target
                    if child.tag.rsplit("}", 1)[-1] == "title"
                )
                if re.search(
                    r"supporting\s+information|supplementary\s+"
                    r"(?:information|materials?|data)",
                    target_title,
                    re.I,
                ) and target not in selected:
                    sections.append(target)
                    selected.add(target)

    # Read actual attachment hrefs from the PMC webpage
    class LinkParser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.links = []

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                href = dict(attrs).get("href")
                if href:
                    self.links.append(href)

    pmc_page_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    pmc_links = None

    def find_attachment_url(filename):
        nonlocal pmc_links

        if pmc_links is None:
            parser = LinkParser()
            parser.feed(read_url(pmc_page_url).decode("utf-8"))
            pmc_links = parser.links

        filename = urllib.parse.unquote(
            urllib.parse.urlsplit(filename).path
        ).rsplit("/", 1)[-1]

        for href in pmc_links:
            link_filename = urllib.parse.unquote(
                urllib.parse.urlsplit(href).path
            ).rsplit("/", 1)[-1]

            if link_filename == filename:
                return urllib.parse.urljoin(pmc_page_url, href)

        raise RuntimeError(
            f"PMC web unable to find supplementary file: {filename}；"
            "The webpage might be blocked by anti-scraping measures, or the attachment link format may require a new rule."
        )

    # Attachment links may be relative filenames inside <media>
    for section in sections:
        source_section = section_name(section)

        for node in section.iter():
            if node.tag.rsplit("}", 1)[-1] != "media":
                continue

            href = next(
                (
                    value for attr, value in node.attrib.items()
                    if attr.rsplit("}", 1)[-1] == "href"
                ),
                "",
            )
            if not href:
                continue

            if re.match(r"https?://", href, re.I):
                add("link", platform_for(href), href, source_section)
            else:
                url = find_attachment_url(href)
                add("link", "PMC", url, source_section)

    # 4. Extract links and raw reads accessions
    accession_pattern = (
        r"\b(?:PRJNA\d+|PRJEB\d+|PRJDB\d+|[SED]R[PRSX]\d+|CRA\d+)\b"
    )

    for section in sections:
        source_section = section_name(section)

        paragraphs = [
            node for node in section.iter()
            if node.tag.rsplit("}", 1)[-1] == "p"
        ]

        for paragraph in paragraphs or [section]:
            urls = []

            for node in paragraph.iter():
                for attr, value in node.attrib.items():
                    if (
                        attr.rsplit("}", 1)[-1] == "href"
                        and re.match(r"https?://", value, re.I)
                    ):
                        urls.append(value)

            content = " ".join(paragraph.itertext())
            urls.extend(re.findall(r"""https?://[^\s<>"']+""", content))
            urls = list(dict.fromkeys(
                url.rstrip(".,;:)") for url in urls
            ))

            has_sra = bool(re.search(
                r"\bSRA\b|\bSequence Read Archive\b", content, re.I
            ))
            raw_urls = []

            for url in urls:
                platform = platform_for(url)
                if (
                    platform == "NCBI"
                    and has_sra
                    and re.search(accession_pattern, urllib.parse.unquote(url), re.I)
                ):
                    add("raw_reads", "NCBI SRA", url, source_section)
                    raw_urls.append(url)
                else:
                    add("link", platform, url, source_section)

            for match in re.finditer(
                accession_pattern, content + " " + " ".join(urls), re.I
            ):
                accession = match.group().upper()

                if any(
                    re.search(
                        r"\b" + re.escape(accession) + r"\b", url, re.I
                    )
                    for url in raw_urls
                ):
                    continue

                if accession.startswith("CRA"):
                    platform = "NGDC GSA"
                    result[:] = [
                        item for item in result
                        if not (
                            item["type"] == "link"
                            and item["platform"] == platform
                            and item["source_section"] == source_section
                            and item["link"].rstrip("/") in {
                                "http://ngdc.cncb.ac.cn/gsa",
                                "https://ngdc.cncb.ac.cn/gsa",
                            }
                        )
                    ]
                elif accession.startswith(("PRJEB", "ER")):
                    platform = "ENA"
                elif accession.startswith(("PRJDB", "DR")):
                    platform = "DDBJ DRA"
                else:
                    # PRJNA is BioProject code, need paragraph shown "SRA" as an evidence
                    if accession.startswith("PRJNA") and not has_sra:
                        continue
                    platform = "NCBI SRA"

                add("raw_reads", platform, accession, source_section)

    return result