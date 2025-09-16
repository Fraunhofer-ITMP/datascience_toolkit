import requests

# Requesting a Token
response = requests.post("https://kggapi.serve.scilifelab.se/request_token", data={"username": "kgg_user1", "password": "gast1@kgg2025"})
token = response.json().get("access_token")
print(token)


# Getting a Disease ID
headers = {"Authorization": f"Bearer {token}"}
response = requests.post("https://kggapi.serve.scilifelab.se/getDiseaseIDs", json={"disease_name": "cancer"}, headers=headers)
print(response.json())


# Creating a KG
headers = {"Authorization": f"Bearer {token}"}
response = requests.post("https://kggapi.serve.scilifelab.se/createKG", json={"disease_id": "MONDO_0004992", "clinical_trial_phase": 3, "protein_threshold": 0.8}, headers=headers)
print(response.json())

