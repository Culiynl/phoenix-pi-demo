import arxiv
import traceback

def test_arxiv():
    try:
        client = arxiv.Client()
        search = arxiv.Search(
            query="cancer",
            max_results=3,
            sort_by=arxiv.SortCriterion.Relevance
        )
        print("Starting search...")
        for r in client.results(search):
            print(f"Title: {r.title}")
            print(f"Authors: {[a.name for a in r.authors]}")
            print(f"Published: {r.published.strftime('%Y-%m-%d')}")
        print("Search complete.")
    except Exception as e:
        traceback.print_exc()

if __name__ == "__main__":
    test_arxiv()
