# OPTIONAL ACQUISITION HELPER. Verify current source terms/permissions before use. Not required for downstream reproduction.
import re
import time
import random
from typing import Optional, Tuple

import pandas as pd
from playwright.sync_api import sync_playwright, TimeoutError as PWTimeoutError
from playwright._impl._errors import TargetClosedError


# ================== 你需要改的参数 ==================
INPUT_XLSX = "中药名单.xlsx"
OUTPUT_XLSX = "中药_性味归经_功能主治_抓取结果.xlsx"
SHEET_NAME = 0
NAME_COL = "中药名"
HEADLESS = True          # 先用 False 调试，成功后改 True
MAX_RETRY = 2
SLEEP_RANGE = (0.6, 1.2)
# 如果你电脑有安装 Google Chrome，强烈建议用它：
USE_SYSTEM_CHROME = True  # True -> channel="chrome"
# ====================================================

MAIN_URL = "https://ydz.chp.org.cn/#/main"


def _clean_text(s: str) -> str:
    s = s.replace("\r\n", "\n")
    s = re.sub(r"[ \t]+\n", "\n", s)
    s = re.sub(r"\n{3,}", "\n\n", s)
    s = re.sub(r"[ \t]{2,}", " ", s)
    return s.strip()


def wait_main_ready(page) -> None:
    """确保 main 页加载完成且搜索框可用。"""
    page.goto(MAIN_URL, wait_until="domcontentloaded", timeout=60000)
    # 等搜索框出现（你截图里 placeholder 就是“请输入关键字”）
    page.get_by_placeholder("请输入关键字").wait_for(state="visible", timeout=30000)


def find_search_input(page):
    """返回可用搜索输入框 locator。"""
    candidates = [
        page.get_by_placeholder("请输入关键字"),
        page.get_by_placeholder(re.compile(r".*(关键字|检索|搜索).*")),
        page.locator("css=input[type='text']").first,
        page.locator("css=input").first,
        page.locator("css=textarea").first,
        page.locator("css=[contenteditable='true']").first,
    ]
    for loc in candidates:
        try:
            if loc.count() == 0:
                continue
            el = loc.first
            el.wait_for(state="visible", timeout=15000)
            el.scroll_into_view_if_needed(timeout=5000)
            return el
        except Exception:
            continue

    page.screenshot(path="debug_no_searchbox.png", full_page=True)
    with open("debug_no_searchbox.html", "w", encoding="utf-8") as f:
        f.write(page.content())
    raise RuntimeError("找不到搜索输入框（已输出 debug_no_searchbox.*）")


def trigger_search(page, box, keyword: str) -> None:
    """输入并触发搜索。"""
    # 先清空再输入
    try:
        box.fill("", timeout=8000)
        box.fill(keyword, timeout=8000)
    except Exception:
        box.click(timeout=8000)
        box.fill("", timeout=8000)
        box.type(keyword, delay=50)

    # 回车触发；不行再点“搜索”按钮
    try:
        box.press("Enter", timeout=3000)
    except Exception:
        btn = page.get_by_role("button", name=re.compile(r"搜索|查询|检索"))
        btn.first.click(timeout=15000)

    # 等结果页/结果列表出现：这里不写死结构，先等页面有变化
    page.wait_for_timeout(1500)


def click_first_result(page) -> bool:
    """
    点击搜索结果第一条进入详情。
    如果站点 DOM 变化，这里是最需要按 debug 图微调的地方。
    """
    candidates = [
        # 常见表格第一行
        page.locator("css=.el-table__body-wrapper tbody tr").first,
        # 常见列表/卡片
        page.locator("css=.result-item").first,
        page.locator("css=.el-list-item").first,
        # 兜底：页面上第一个看起来像“条目链接”的 a
        page.locator("xpath=(//a[contains(@href,'#/')])[1]"),
    ]
    for loc in candidates:
        try:
            if loc.count() == 0:
                continue
            loc.scroll_into_view_if_needed(timeout=5000)
            loc.click(timeout=15000)
            page.wait_for_timeout(1500)
            return True
        except Exception:
            continue
    return False


def extract_section(page, heading: str):
    """
    适配 2020版药典页面结构：
    标题为： 【性味与归经】
    内容为： 下一个或后续连续的 <p> 段落
    """

    # 1. 找到标题段落
    title = page.locator(f"xpath=//p[normalize-space()='{heading}' or normalize-space()='【{heading}】']").first

    if title.count() == 0:
        return None

    # 2. 获取标题后的所有兄弟段落
    paragraphs = title.locator("xpath=following-sibling::p")

    texts = []

    for i in range(paragraphs.count()):
        p = paragraphs.nth(i)
        text = p.inner_text().strip()

        # 遇到下一个章节标题就停止
        if text.startswith("【") and text.endswith("】"):
            break

        if text:
            texts.append(text)

    if not texts:
        return None

    return "\n".join(texts)


def is_detail_page(page) -> bool:
    # 详情页一定会同时出现这些“章节标题”中的至少一个
    # 且页面不应只是一堆“数据库搜索/检索”之类的导航词
    has_section = (
        page.locator("text=性味与归经").count() > 0
        or page.locator("text=功能与主治").count() > 0
        or page.locator("text=性味归经").count() > 0
        or page.locator("text=功能主治").count() > 0
    )
    if not has_section:
        return False

    # 进一步排除：如果页面明显是结果/检索页（常见关键词）
    bad_words = ["数据库搜索", "检索结果", "搜索结果", "目录", "返回", "筛选"]
    body_text = page.locator("body").inner_text(timeout=5000)
    hit_bad = any(w in body_text for w in bad_words)

    # 命中 bad_words 并不一定是错（详情页也可能有“返回”），
    # 但如果命中很多，就判为非详情页
    bad_count = sum(1 for w in bad_words if w in body_text)
    return bad_count <= 1

def open_detail_from_results(page, herb_name: str):
    # 返回：详情页 page 或 None
    herb_name = str(herb_name).strip()
    if not herb_name:
        return None

    candidates = [
        page.locator(f"xpath=//*[self::a or self::div or self::span][normalize-space()='{herb_name}']").first,
        page.locator(f"xpath=(//a[contains(normalize-space(), '{herb_name}')])[1]").first,
        page.locator("css=.el-table__body-wrapper tbody tr").first,
        page.locator("xpath=(//a[contains(@href,'#/')])[1]").first,
    ]

    for loc in candidates:
        try:
            if loc.count() == 0:
                continue
            loc.scroll_into_view_if_needed(timeout=5000)

            # popup 新页
            try:
                with page.expect_popup(timeout=2000) as pop:
                    loc.click(timeout=8000)
                new_page = pop.value
                new_page.wait_for_load_state("domcontentloaded", timeout=30000)
                page = new_page
            except PWTimeoutError:
                loc.click(timeout=8000)

            for _ in range(30):
                page.wait_for_timeout(300)
                if is_detail_page(page):
                    return page
        except Exception:
            continue

    return None

def fetch_one(page, name: str):
    name = str(name).strip()
    if not name:
        return None, None, "空名称"

    page.goto("https://ydz.chp.org.cn/#/main", wait_until="domcontentloaded", timeout=60000)

    box = page.get_by_placeholder("请输入关键字")
    box.wait_for(state="visible", timeout=20000)
    box.fill("")
    box.fill(name)
    page.get_by_text("搜索", exact=True).click()
    page.wait_for_timeout(1200)

    # 如果搜索后不是详情页，就从结果页点进去（兼容新标签页）
    if not is_detail_page(page):
        detail_page = open_detail_from_results(page, name)  # 用“返回page版本”
        if detail_page is None:
            page.screenshot(path=f"debug_{name}_still_results.png", full_page=True)
            with open(f"debug_{name}_still_results.html", "w", encoding="utf-8") as f:
                f.write(page.content())
            return None, None, "停留在结果页，无法进入详情"
        page = detail_page  # 关键：切换到真正的详情页

    # 再确认一次
    if not is_detail_page(page):
        page.screenshot(path=f"debug_{name}_not_detail.png", full_page=True)
        return None, None, "未进入详情页"

    xw = extract_section(page, "性味与归经") or extract_section(page, "性味归经")
    gn = extract_section(page, "功能与主治") or extract_section(page, "功能主治")

    # 脏数据过滤
    bad = {"数据库搜索", "搜索", "检索结果", "搜索结果", "目录"}
    def ok_text(t: str) -> bool:
        if not t:
            return False
        t = t.strip()
        if t in bad:
            return False
        if len(t) < 6:  # 太短一般不是正文
            return False
        return True

    xw = xw if ok_text(xw) else None
    gn = gn if ok_text(gn) else None

    if not xw and not gn:
        page.screenshot(path=f"debug_{name}_nofield.png", full_page=True)
        with open(f"debug_{name}_nofield.html", "w", encoding="utf-8") as f:
            f.write(page.content())
        return None, None, "进详情但字段未抓到或为脏数据"

    return xw, gn, "OK"


def main():
    df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)
    if NAME_COL not in df.columns:
        raise ValueError(f"Excel 中找不到列：{NAME_COL}。现有列：{list(df.columns)}")

    for col in ["性味与归经", "功能与主治", "抓取状态"]:
        if col not in df.columns:
            df[col] = ""

    with sync_playwright() as p:
        launch_kwargs = dict(
            headless=HEADLESS,
            args=[
                "--start-maximized",
                "--ignore-certificate-errors",
                "--disable-blink-features=AutomationControlled",
            ],
        )
        if USE_SYSTEM_CHROME:
            launch_kwargs["channel"] = "chrome"

        browser = p.chromium.launch(**launch_kwargs)

        context = browser.new_context(
            ignore_https_errors=True,
            viewport=None,
            locale="zh-CN",
        )
        page = context.new_page()

        try:
            total = len(df)
            for i, row in df.iterrows():
                name = row[NAME_COL]
                if str(row.get("抓取状态", "")).strip() == "OK":
                    continue

                last_err = None
                for attempt in range(1, MAX_RETRY + 2):
                    try:
                        xw, gn, st = fetch_one(page, name)
                        df.at[i, "性味与归经"] = xw or ""
                        df.at[i, "功能与主治"] = gn or ""
                        df.at[i, "抓取状态"] = st
                        print(f"[{i+1}/{total}] {name} -> {st}")
                        break
                    except TargetClosedError:
                        # 页面被关闭就新开一个
                        page = context.new_page()
                        last_err = "TargetClosedError"
                        continue
                    except Exception as e:
                        last_err = e
                        print(f"[{i+1}/{total}] {name} -> 第{attempt}次异常：{type(e).__name__}: {e}，重试中…")
                        page.wait_for_timeout(800)

                if str(df.at[i, "抓取状态"]).strip() == "":
                    df.at[i, "抓取状态"] = f"异常: {type(last_err).__name__}: {last_err}" if last_err else "未知异常"

                time.sleep(random.uniform(*SLEEP_RANGE))

        except KeyboardInterrupt:
            print("检测到手动中断（Ctrl+C），保存已抓取结果…")

        finally:
            df.to_excel(OUTPUT_XLSX, index=False)
            print(f"\n已保存：{OUTPUT_XLSX}")
            context.close()
            browser.close()


if __name__ == "__main__":
    main()

#如果你发现功能匹配“太严格/太宽松”，优先调这两个参数就行：

#RATIO_THR（比如 0.70 更严格，0.60 更宽松）

#BIGRAM_THR（比如 0.65 更严格，0.50 更宽松）