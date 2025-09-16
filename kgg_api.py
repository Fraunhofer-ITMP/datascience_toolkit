"""FastAPI for KGG (Knowledge Graph Generation). This FastAPI app corresponds to pages/2_KGG.py in the streamlit app."""

import datetime
import os
import pandas as pd
import logging
import uvicorn
from dotenv import load_dotenv
import pybel
from fastapi import FastAPI, HTTPException, Depends, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from jose import jwt, JWTError
import datetime
from api_utils import kgg_apiutils

load_dotenv()
app = FastAPI()

origins = [
    "http://localhost:8080",
    "0.0.0.0:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class KGGCreation(BaseModel):
    disease_id: str = "MONDO_0004976"
    clinical_trial_phase: int = 3
    protein_threshold: float = 0.8
    created_kg: dict = None

class DiseaseID(BaseModel):
    disease_name: str = "cancer"
    disease_ids: list = None

class Token(BaseModel):
    access_token: str
    token_type: str


####### SECURITY RELATED FUNCTIONS ########################
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")

def validate_user(username: str, password: str):
    """User validation using environment variables"""
    valid_users = {
        os.getenv("KGG_API_USER1_USERNAME"): os.getenv("KGG_API_USER1_PASSWORD"),
    }
    return valid_users.get(username) == password

def get_current_user(token: str = Depends(oauth2_scheme)):
    try:
        payload = jwt.decode(token, os.getenv("KGG_API_SECRET"), algorithms=[os.getenv("KGG_SECRET_ALGORITHM")])
        return payload["sub"]
    except JWTError:
        raise HTTPException(status_code=401, detail="Invalid token")
###########################################################

################### LOGGING ###############################
def logging_setup(request):
    request_ip_address = request.client.host
    request_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    if not os.path.exists('kgg_logs'):
        os.makedirs('kgg_logs')
    log_filename = f"kgg_logs/{request_ip_address}_{request_time.replace(':', '-')}.log"
    
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
##########################################################


@app.post("/request_token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    if not validate_user(form_data.username, form_data.password):
        raise HTTPException(
            status_code=401,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
#    access_token_expires = timedelta(days=30)          # for now not expiring anything ;)
    access_token = jwt.encode(
        {
            "sub": form_data.username,
#            "exp": datetime.utcnow() + access_token_expires
        },
        os.getenv("KGG_API_SECRET"),
        algorithm=os.getenv("KGG_SECRET_ALGORITHM")
    )
    
    return {"access_token": access_token, "token_type": "bearer"}


@app.get("/", response_class=HTMLResponse)
async def read_root():
    """Root endpoint."""
    html_content = """
    <html>
        <head>
            <title>Knowledge Generator API</title>
        </head>
        <body>
        <h1>Endpoints:</h1>
        1. <a href="/createKG">/createKG</a>: <br>Generates a Knowledge Graph based on disease ID, clinical trial phase, and protein threshold.<br>

                1.1 Input Format:
                <pre>
                {
                    "disease_id": "MONDO_0004976",
                    "clinical_trial_phase": 3,
                    "protein_threshold": 0.8
                }
                </pre>
                1.2 Output Format:
                <pre>
                {
                    "disease_id": "MONDO_0004976",
                    "clinical_trial_phase": 3,
                    "protein_threshold": 0.8,
                    "created_kg": { ... }  # JGIF JSON format of the generated Knowledge Graph. Please look at PyBel documentation (pybel.from_jgif_json) for more details on how to parse this format.
                }
                </pre>

        2. <a href="/getDiseaseIDs">/getDiseaseIDs</a>: <br>Returns available disease IDs from a given disease name.<br>
        
                2.1 Input Format:
                <pre>
                {
                    "disease_name": "cancer"
                }
                </pre>
                2.2 Output Format:
                <pre>
                {
                    "disease_ids": A list of dictionary containing disease IDs and other information.
                }
                </pre>

        <h2>Authentication</h2>
        <p> First you need to get a token by calling the <a href="/token">/request_token</a> endpoint with your username and password. <br> 
        (To get the username and password, please contact Fraunhofer ITMP).
        Once you have a username and a password, you can generate token for APIs in following way in python:</p>
        <pre>
        import requests

        response = requests.post("http://10.164.197.141:8080/request_token", data={"username": "YOUR_USERNAME", "password": "YOUR_PASSWORD"})
        token = response.json().get("access_token")
        </pre>
        <p>Now you can use this token in the Authorization header as a Bearer token for all other endpoints.</p>

        <p>Example of using the token in python:</p>
        <pre>
        headers = {"Authorization": f"Bearer {token}"}
        response = requests.post("http://10.164.197.141:8080/createKG", json={"disease_id": "MONDO_0004976", "clinical_trial_phase": 3, "protein_threshold": 0.8}, headers=headers)
        print(response.json())
        </pre>
        Please follow this same format for all other endpoints that might need authentication!

    </body>
    </html>

    """
    return HTMLResponse(content=html_content, status_code=200)

@app.post("/getDiseaseIDs",response_model=DiseaseID)
async def get_disease_ids(disease_model:DiseaseID, request: Request, current_user: str = Depends(get_current_user)):
    """This function returns available disease IDs using searchDisease function from kgg_apiutils.py"""
    logging_setup(request)    
    logging.info(f"Current user: {current_user}")
    try:
        disease_name = disease_model.disease_name.strip()
        if not disease_name:
            logging.error("Disease name is required.")
            raise HTTPException(status_code=400, detail="Disease name is required.")
        disease_ids = kgg_apiutils.searchDisease(keyword=disease_name,logger=logging.getLogger())
        logging.info(f"Total disease IDs found: {len(disease_ids)} for keyword '{disease_name}'")
        if disease_ids.empty:
            logging.warning("No disease IDs found.")
            raise HTTPException(status_code=404, detail="No disease IDs found.")
        
        disease_ids_list = disease_ids.to_dict('records')       # Converting dataframe to a list of dictionaries so that it can be returned as JSON
        logging.info(f"Disease IDs: {disease_ids_list}")
        return {"disease_ids": disease_ids_list}    
    except Exception as e:
        logging.error(f"An error occurred while fetching disease IDs: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/createKG", response_model=KGGCreation)
async def generate_kg(kgg_model: KGGCreation, request: Request, current_user: str = Depends(get_current_user)):
    """Generate a Knowledge Graph based on the provided parameters."""
    logging_setup(request)
    
    logging.info(f"Received parameters: {kgg_model.dict()}")
    logging.info(f"Current user: {current_user}")

    try:
        disease_id = kgg_model.disease_id.strip()
        clinical_trial_phase = kgg_model.clinical_trial_phase
        protein_threshold = kgg_model.protein_threshold
        if not disease_id:
            logging.error("Disease name is required.")
            raise HTTPException(status_code=400, detail="Disease name is required.")
        if not clinical_trial_phase:
            logging.info("Clinical trial phase not provided, defaulting to 1.")
            clinical_trial_phase = 1
        if not protein_threshold:
            logging.info("Protein threshold not provided, defaulting to 0.0.")
            protein_threshold = 0.0

        created_kg = kgg_apiutils.createKG(disease_id=disease_id,
                                        clinical_trial_phase=clinical_trial_phase,
                                        protein_threshold=protein_threshold,
                                        logger=logging.getLogger())

        logging.info("Knowledge Graph generation completed successfully.")
        kgg_model.created_kg = pybel.to_jgif_jsons(created_kg) if created_kg else None
        return JSONResponse(content=kgg_model.dict(), status_code=200)
    except HTTPException as http_exception:
        logging.error(f"HTTP Exception: {http_exception.detail}")
        raise http_exception
    
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))
    
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080, log_level="info")