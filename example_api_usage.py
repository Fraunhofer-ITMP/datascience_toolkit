import requests

# Requesting a Token
response = requests.post("http://localhost:8080/request_token", data={"username": "kgg_user1", "password": "gast1@kgg2025"})
token = response.json().get("access_token")
print(token)


# Getting a Disease ID
headers = {"Authorization": f"Bearer {token}"}
response = requests.post("http://localhost:8080/getDiseaseIDs", json={"disease_name": "covid"}, headers=headers)
print(response.json())


# Creating a KG
headers = {"Authorization": f"Bearer {token}"}
response = requests.post("http://localhost:8080/createKG", json={"disease_id": "MONDO_0100233", "clinical_trial_phase": 3, "protein_threshold": 0.8}, headers=headers)
print(response.json())

