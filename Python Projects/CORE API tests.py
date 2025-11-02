import requests

import json

from docling.document_converter import DocumentConverter

import pymupdf4llm
"""
def search_outputs_by_keywords(keywords:str):
    url = "https://api.core.ac.uk/v3/search/outputs"
    headers = {"Authorization": f"Bearer {"zypYaFWSBqE8Ngef5AjI9li36PvnUcdo"}"}
    
    params = {"q": f"fullText:{keywords} AND _exists_:downloadUrl AND language.code:en ", "limit": 1}
    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            print(response.content)
            return response.json()
        else:
            print(f"Error {response.status_code}: {response.text}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return None

data = search_outputs_by_keywords("Fischer Tropsch")
    
with open('data.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
    
results_list = data.get("results")

download_urls = []

for results in results_list:
    download_urls.append(results["downloadUrl"][0])

print(download_urls[0])
source = download_urls[0]
"""
#pymupdf4llm example
#md_text = pymupdf4llm.to_markdown(source)
#print(md_text)

#docling example
converter = DocumentConverter()
result = converter.convert(f"C:/Users/wille/Downloads/A COMPREHENSIVE REVIEW OF FISCHER-TROPSCH SYNTHESIS PROCESSES TO.pdf")
markdown = (result.document.export_to_markdown())

output_markdown_file_docling = 'output_docling.md'
with open(output_markdown_file_docling, 'w', encoding='utf-8') as f:
    f.write(markdown)
print(f"Markdown from docling saved to {output_markdown_file_docling}")
    
    