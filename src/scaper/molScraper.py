from bs4 import BeautifulSoup
import requests
import re

# Parse the HTML
url = """https://www.genome.jp/dbget-bin/www_bget?cpd:C00819"""  # Your HTML here
response = requests.get(url)
soup = BeautifulSoup(response.content, 'html.parser')

# Find the button with the text "Mol file"
button = soup.find('button', text='Mol file')

# Extract the URL from the onclick attribute
onclick = button['onclick']
url = re.search(r"'(.*?)'", onclick).group(1)  # Extract the URL between single quotes

# Add the base URL if the URL is relative
if not url.startswith('http'):
    url = 'https://www.genome.jp/dbget-bin/www_bget?cpd:C00819/' + url

# Send a GET request to the URL
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    # Write the content to a file
    with open('mol_file.mol', 'wb') as file:
        file.write(response.content)