# OPTIONAL ACQUISITION HELPER. Verify current source terms/permissions before use. Not required for downstream reproduction.
import time
from typing import Dict, List

import pandas as pd
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm


BASE = "https://bidd.group/TCMID"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/122.0 Safari/537.36"
    )
}


def get_soup(url: str, session: requests.Session, timeout: int = 30) -> BeautifulSoup:
    r = session.get(url, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    return BeautifulSoup(r.text, "lxml")


def parse_targets_from_component_page(session: requests.Session, tcmh_id: str) -> List[Dict]:
    """
    从 herb.php?herb=TCMHxxxx 页面解析靶点表
    返回字段：
      TCMH_ID, Target_ID, Gene_Symbol, Target_Name, Target_Class, Uniprot_ID, Source_URL
    """
    url = f"{BASE}/herb.php?herb={tcmh_id}"
    soup = get_soup(url, session)

    # 找包含 Targeted Human Proteins 的区域附近表格
    anchor = None
    for tag in soup.find_all(["h1", "h2", "h3", "h4", "h5", "div", "p"]):
        txt = tag.get_text(" ", strip=True)
        if "Targeted Human Proteins" in txt:
            anchor = tag
            break

    candidate_tables = []
    if anchor:
        candidate_tables = list(anchor.find_all_next("table"))
    else:
        candidate_tables = soup.find_all("table")

    results = []

    for table in candidate_tables:
        headers = [th.get_text(" ", strip=True) for th in table.find_all("th")]
        if not headers:
            continue

        header_join = " | ".join(headers).lower()
        if ("target id" in header_join) and ("target name" in header_join):
            for tr in table.find_all("tr"):
                tds = [td.get_text(" ", strip=True) for td in tr.find_all("td")]
                if len(tds) < 2:
                    continue

                row = {headers[i]: tds[i] if i < len(tds) else "" for i in range(len(headers))}

                def pick(*keys):
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

            if results:
                break

    return results


def normalize_str(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def main(
    input_xlsx: str = "TCMID_靶点抓取结果.xlsx",
    output_xlsx: str = "TCMID_靶点抓取结果_补充后.xlsx",
    sleep_s: float = 0.3
):
    # 1) 读取两个 sheet
    df_targets = pd.read_excel(input_xlsx, sheet_name="targets_long")
    df_map = pd.read_excel(input_xlsx, sheet_name="name_to_tcmh_mapping")

    # 统一清洗
    df_targets["Herb_ChineseName"] = df_targets["Herb_ChineseName"].map(normalize_str)
    df_targets["TCMH_ID"] = df_targets["TCMH_ID"].map(normalize_str)

    df_map["Herb_ChineseName"] = df_map["Herb_ChineseName"].map(normalize_str)
    df_map["Chosen_TCMH_ID"] = df_map["Chosen_TCMH_ID"].map(normalize_str)

    # 2) 建立 targets_long 中已经存在的 “药物名 -> TCMH_ID集合”
    existed_map = (
        df_targets.groupby("Herb_ChineseName")["TCMH_ID"]
        .apply(lambda s: set(i for i in s if i))
        .to_dict()
    )

    # 3) 从 sheet2 找出需要补充的记录
    to_supplement = []
    for _, row in df_map.iterrows():
        herb = row["Herb_ChineseName"]
        chosen_id = row["Chosen_TCMH_ID"]

        if not herb or not chosen_id:
            continue

        existed_ids = existed_map.get(herb, set())

        # 只补 sheet1 里还没有的药物/ID 组合
        if chosen_id not in existed_ids:
            to_supplement.append((herb, chosen_id))

    print(f"需补充的药物数量: {len(to_supplement)}")

    # 4) 抓取补充数据
    new_rows = []
    failed_rows = []

    with requests.Session() as session:
        for herb, tcmh_id in tqdm(to_supplement, desc="补充抓取中"):
            try:
                targets = parse_targets_from_component_page(session, tcmh_id)

                if targets:
                    for t in targets:
                        t["Herb_ChineseName"] = herb
                        new_rows.append(t)
                else:
                    # 没有靶点也保留一行，避免下次重复抓
                    new_rows.append({
                        "TCMH_ID": tcmh_id,
                        "Target_ID": "",
                        "Gene_Symbol": "",
                        "Target_Name": "",
                        "Target_Class": "",
                        "Uniprot_ID": "",
                        "Source_URL": f"{BASE}/herb.php?herb={tcmh_id}",
                        "Herb_ChineseName": herb,
                    })

                time.sleep(sleep_s)

            except Exception as e:
                failed_rows.append({
                    "Herb_ChineseName": herb,
                    "TCMH_ID": tcmh_id,
                    "Error": str(e)
                })
                print(f"[WARN] {herb} - {tcmh_id} 抓取失败: {e}")

    # 5) 合并旧数据和新数据
    df_new = pd.DataFrame(new_rows)

    if not df_new.empty:
        df_all = pd.concat([df_targets, df_new], ignore_index=True)

        # 去重：按药物名 + TCMH_ID + Target_ID 去重
        # 若 Target_ID 为空，也能保留一条“空靶点占位记录”
        df_all["_dedup_key"] = (
            df_all["Herb_ChineseName"].map(normalize_str) + "||" +
            df_all["TCMH_ID"].map(normalize_str) + "||" +
            df_all["Target_ID"].map(normalize_str)
        )
        df_all = df_all.drop_duplicates(subset=["_dedup_key"], keep="first").drop(columns=["_dedup_key"])
    else:
        df_all = df_targets.copy()

    df_failed = pd.DataFrame(failed_rows)

    # 6) 导出
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        df_all.to_excel(writer, index=False, sheet_name="targets_long")
        df_map.to_excel(writer, index=False, sheet_name="name_to_tcmh_mapping")
        if not df_failed.empty:
            df_failed.to_excel(writer, index=False, sheet_name="failed_log")

    print(f"补充完成，输出文件：{output_xlsx}")
    print(f"新增记录数：{len(df_new)}")
    print(f"失败记录数：{len(df_failed)}")


if __name__ == "__main__":
    main()