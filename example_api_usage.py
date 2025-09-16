import requests

# Requesting a Token
response = requests.post("http://10.164.197.141:8080/request_token", data={"username": "kgg_user1", "password": "gast1@kgg2025"})
token = response.json().get("access_token")
print(token)

# Making an API call
headers = {"Authorization": f"Bearer {token}"}
response = requests.post("http://10.164.197.141:8080/createKG", json={"disease_id": "MONDO_0004976", "clinical_trial_phase": 3, "protein_threshold": 0.8}, headers=headers)
print(response.json())

