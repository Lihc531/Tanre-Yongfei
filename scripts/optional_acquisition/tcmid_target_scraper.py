# OPTIONAL ACQUISITION HELPER. Verify current source terms/permissions before use. Not required for downstream reproduction.
import re
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm


BASE = "https://bidd.group/TCMID"
BROWSE_URL = f"{BASE}/browse.php"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/122.0 Safari/537.36"
}


def get_soup(url: str, session: requests.Session, timeout: int = 30) -> BeautifulSoup:
    """Fetch url and return BeautifulSoup."""
    r = session.get(url, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    return BeautifulSoup(r.text, "lxml")


def extract_funcclass_links(session: requests.Session) -> List[str]:
    """
    From browse.php, extract all links like results.php?browsefuncclass=...
    These pages contain component tables with (Component ID, Name-Chinese, ...).
    """
    soup = get_soup(BROWSE_URL, session)
    links = []
    for a in soup.select("a[href]"):
        href = a.get("href", "")
        if "results.php?browsefuncclass=" in href:
            if href.startswith("http"):
                links.append(href)
            else:
                links.append(f"{BASE}/{href.lstrip('/')}")
    # 去重（保持顺序）
    seen = set()
    uniq = []
    for x in links:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return uniq


def parse_component_table(page_soup: BeautifulSoup) -> List[Tuple[str, str]]:
    """
    Parse the 'TCM Component' table in results.php pages.
    Return list of (TCMH_id, chinese_name).
    """
    pairs = []

    # 经验上这些 results.php 页面里会有一个组件表格，且行内含 TCMHxxxx
    for row in page_soup.select("tr"):
        cols = [c.get_text(" ", strip=True) for c in row.select("td")]
        if not cols:
            continue

        # 常见列顺序：Component ID | Name-PinYin | Name-Chinese | Name-English | Barcode
        # 我们抓：TCMHxxxx 和中文名
        tcmh = None
        zh = None

        for c in cols:
            m = re.search(r"\bTCMH\d+\b", c)
            if m:
                tcmh = m.group(0)
                break

        if tcmh:
            # 尝试从 cols 里找“中文名”那一列：一般是第3列（索引2）
            if len(cols) >= 3:
                candidate = cols[2]
                # 中文名列往往含中文字符
                if re.search(r"[\u4e00-\u9fff]", candidate):
                    zh = candidate

            # 如果没取到，用整行里第一个“含中文”的列顶上
            if not zh:
                for c in cols:
                    if re.search(r"[\u4e00-\u9fff]", c):
                        zh = c
                        break

        if tcmh and zh:
            pairs.append((tcmh, zh))

    return pairs


def build_cn_to_tcmh_map(session: requests.Session, sleep_s: float = 0.2) -> Dict[str, List[str]]:
    """
    Crawl all functional class pages and build mapping:
      ChineseName -> [TCMH... , ...]
    """
    func_links = extract_funcclass_links(session)
    if not func_links:
        raise RuntimeError("未能从 browse.php 抓到 browsefuncclass 链接（可能网络/站点结构变化）。")

    cn2ids: Dict[str, set] = {}
    for url in tqdm(func_links, desc="Crawling functional classes"):
        try:
            soup = get_soup(url, session)
            pairs = parse_component_table(soup)
            for tcmh, zh in pairs:
                cn2ids.setdefault(zh, set()).add(tcmh)
        except Exception as e:
            # 不中断全局任务，记录一下
            print(f"[WARN] Failed: {url} -> {e}")

        time.sleep(sleep_s)

    # set -> sorted list
    cn2ids_list: Dict[str, List[str]] = {k: sorted(list(v)) for k, v in cn2ids.items()}
    return cn2ids_list


def choose_best_id(ids: List[str]) -> str:
    """
    If multiple IDs for same Chinese name, choose a stable one:
    - Prefer the smallest numeric id (often the 'main' entry).
    """
    def key_fn(x: str) -> int:
        m = re.search(r"TCMH(\d+)", x)
        return int(m.group(1)) if m else 10**18

    return sorted(ids, key=key_fn)[0]


def parse_targets_from_component_page(session: requests.Session, tcmh_id: str) -> List[Dict]:
    """
    Fetch herb.php?herb=TCMHxxxx and parse target table.
    Output list of dict rows with:
      Target ID, Gene Symbol, Target Name, Target Class, Uniprot ID
    """
    url = f"{BASE}/herb.php?herb={tcmh_id}"
    soup = get_soup(url, session)

    # 找到包含 “Targeted Human Proteins” 的标题附近的表格
    # 页面结构可能变化，因此采用“文本锚点 + 向后找表格”的方式
    anchor = None
    for tag in soup.find_all(["h1", "h2", "h3", "h4", "h5", "p", "div"]):
        txt = tag.get_text(" ", strip=True)
        if "Targeted Human Proteins" in txt:
            anchor = tag
            break

    tables = soup.find_all("table")
    if not tables:
        return []

    # 如果找到 anchor：从 anchor 后面最近的 table 开始尝试解析
    candidate_tables = []
    if anchor:
        # 取 anchor 之后的所有 table
        for t in anchor.find_all_next("table"):
            candidate_tables.append(t)
    else:
        candidate_tables = tables  # 兜底

    results = []
    for t in candidate_tables:
        # 读表头
        headers = [th.get_text(" ", strip=True) for th in t.find_all("th")]
        if not headers:
            continue

        # 目标表常见表头包含：Target ID / Gene Symbol / Target Name / Target Class / Uniprot
        header_join = " | ".join(headers).lower()
        if ("target id" in header_join) and ("target name" in header_join):
            # 解析每行
            for tr in t.find_all("tr"):
                tds = [td.get_text(" ", strip=True) for td in tr.find_all("td")]
                if len(tds) < 2:
                    continue

                row = {headers[i]: tds[i] if i < len(tds) else "" for i in range(len(headers))}
                # 归一化字段名（不同页面可能大小写/空格不同）
                def pick(*keys) -> str:
                    for k in keys:
                        for real_k in row.keys():
                            if real_k.strip().lower() == k.strip().lower():
                                return row.get(real_k, "")
                    return ""

                target_id = pick("Target ID")
                target_name = pick("Target Name")
                if not target_id and not target_name:
                    continue

                results.append({
                    "TCMH_ID": tcmh_id,
                    "Target_ID": target_id,
                    "Gene_Symbol": pick("Gene Symbol", "Gene"),
                    "Target_Name": target_name,
                    "Target_Class": pick("Target Class", "Class"),
                    "Uniprot_ID": pick("Uniprot ID", "Uniprot", "UniProt ID"),
                    "Source_URL": url,
                })

            # 找到一个合格的目标表就够了（避免重复解析）
            if results:
                break

    return results


def main(
    input_xlsx: str = "筛选后的传播模型候选中药.xlsx",#读取文件
    herb_col: str = "中药名",
    output_xlsx: str = "TCMID_靶点抓取结果.xlsx",
    sleep_s: float = 0.2
):
    with requests.Session() as session:
        # 1) 读输入
        df = pd.read_excel(input_xlsx)
        if herb_col not in df.columns:
            raise ValueError(f"输入表格缺少列：{herb_col}；当前列有：{list(df.columns)}")

        herbs = (
            df[herb_col]
            .astype(str)
            .str.strip()
            .replace({"nan": None, "None": None, "": None})
            .dropna()
            .unique()
            .tolist()
        )

        # 2) 建映射：中文名 -> TCMH ids
        print("Step 1/3: Building ChineseName -> TCMH mapping from browse pages ...")
        cn2ids = build_cn_to_tcmh_map(session, sleep_s=sleep_s)

        # 3) 逐个药名匹配 TCMH 并抓 targets
        print("Step 2/3: Fetching targets for each herb ...")
        all_rows = []
        mapping_rows = []

        for zh_name in tqdm(herbs, desc="Herbs"):
            ids = cn2ids.get(zh_name, [])
            chosen_id = choose_best_id(ids) if ids else None
            mapping_rows.append({
                "Herb_ChineseName": zh_name,
                "Matched_TCMH_IDs": ",".join(ids) if ids else "",
                "Chosen_TCMH_ID": chosen_id or "",
            })

            if not chosen_id:
                continue

            try:
                targets = parse_targets_from_component_page(session, chosen_id)
                if not targets:
                    # 也输出一行空靶点（方便你知道抓到了但没靶点）
                    all_rows.append({
                        "Herb_ChineseName": zh_name,
                        "TCMH_ID": chosen_id,
                        "Target_ID": "",
                        "Gene_Symbol": "",
                        "Target_Name": "",
                        "Target_Class": "",
                        "Uniprot_ID": "",
                        "Source_URL": f"{BASE}/herb.php?herb={chosen_id}",
                    })
                else:
                    for t in targets:
                        t["Herb_ChineseName"] = zh_name
                        all_rows.append(t)

            except Exception as e:
                print(f"[WARN] Failed to parse targets for {zh_name} ({chosen_id}): {e}")

            time.sleep(sleep_s)

        # 4) 输出
        print("Step 3/3: Writing output Excel ...")
        out_targets = pd.DataFrame(all_rows)
        out_map = pd.DataFrame(mapping_rows)

        with pd.ExcelWriter(output_xlsx, engine="openpyxl") as w:
            out_targets.to_excel(w, index=False, sheet_name="targets_long")
            out_map.to_excel(w, index=False, sheet_name="name_to_tcmh_mapping")

        print(f"Done! Output: {output_xlsx}")


if __name__ == "__main__":
    main()