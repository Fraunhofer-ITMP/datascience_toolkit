"""FastAPI for KGG (Knowledge Graph Generation). This FastAPI app corresponds to pages/2_KGG.py in the streamlit app."""

import datetime
import os

import logging
import uvicorn
from dotenv import load_dotenv

from fastapi import FastAPI, HTTPException, Depends, status, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.security import OAuth2PasswordBearer
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from jose import jwt, JWTError

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

class KGG(BaseModel):
    disease: str = "COVID-19"
    clinical_trial_phase: int = 1
    protein_threshold: float = 0.0



####### SECURITY RELATED FUNCTIONS ########################
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")
def get_current_user(token: str = Depends(oauth2_scheme)):
    try:
        payload = jwt.decode(token, os.getenv("KGG_API_SECRET"), algorithms=[os.getenv("KGG_SECRET_ALGORITHM")])
        return payload["sub"]
    except JWTError:
        raise HTTPException(status_code=401, detail="Invalid token")
###########################################################

################### LOGGING ###############################
def logging_setup(Request):
    request_ip_address = Request.client.host
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


############### KGG RELATED FUNCTIONS ####################
def process_disease(disease_id: str, clinical_trial_phase: int, protein_threshold: float):
    """Performs everything that 2_KGG.py does."""
    disease_df = kgg_utils.searchDisease(disease.strip())
    logging.info(f"Searching for disease: {disease}")
    if disease_df.empty:
        logging.error(f"Disease Dataframe is not created because disease is not found for: {disease}")
        raise HTTPException(status_code=404, detail=f"No data found for disease: {disease}")
    disease_id = disease_df.iloc[0]["id"]
    logging.info(f"Disease ID found: {disease_id}")
#    viral_proteins = kgg_utils.GetViralProteins(disease)
    # logging.info(f"{viral_proteins}")
    # if viral_proteins:
    #     logging.info(f"Viral proteins found: {len(viral_proteins)}")
    # else:
    #     logging.info(f"Viral Protein not found.")
    drugs_df, dis2prot_df, dis2snp_df = kgg_utils.createInitialKG(_ct_phase=clinical_trial_phase)
    if drugs_df or dis2prot_df or dis2snp_df:
        logging.info(f"Viral proteins found: {len(viral_proteins)}")

    filtered_dis2prot_df = dis2prot_df[dis2prot_df["Score"] >= protein_threshold]
    json_output = {
        "disease": disease,
#        "viral_proteins": viral_proteins,
        "drugs_df": drugs_df,
        "dis2prot_df": filtered_dis2prot_df,
        "dis2snp_df": dis2snp_df
    }
    logging.info(f"Initial Knowledge Graph created for disease: {disease}")
    return json_output


@app.get("/", response_class=HTMLResponse)
async def read_root():
    """Root endpoint."""
    return HTMLResponse(content="<h1>Knowledge Generator API</h1> Please follow the following instructions to use the API: "
    "<br> 1. Use the /generateKG endpoint to generate a Knowledge Graph.<br>" 
    "1.1 Parameters required {'disease': 'YOUR_DISEASE_NAME'} <br> ", status_code=200)

@app.post("/generateKG", response_model=KGG)
async def generate_kg(kgg: KGG, request: Request):
    """Generate a Knowledge Graph based on the provided parameters."""
    logging_setup(request)
    
    logging.info(f"Received parameters: {kgg.dict()}")
    
    # Checking security
    current_user = get_current_user()
    logging.info(f"Current user: {current_user}")

    try:
        disease = kgg.disease.strip()
        clinical_trial_phase = kgg.clinical_trial_phase
        protein_threshold = kgg.protein_threshold
        if not disease:
            logging.error("Disease name is required.")
            raise HTTPException(status_code=400, detail="Disease name is required.")
        if not clinical_trial_phase:
            logging.info("Clinical trial phase not provided, defaulting to 1.")
            clinical_trial_phase = 1
        if not protein_threshold:
            logging.info("Protein threshold not provided, defaulting to 0.0.")
            protein_threshold = 0.0


        process_disease(disease, clinical_trial_phase, protein_threshold)
        logging.info("Knowledge Graph generation completed successfully.")
        return JSONResponse(content=kgg.dict(), status_code=200)
    except HTTPException as http_exception:
        logging.error(f"HTTP Exception: {http_exception.detail}")
        raise http_exception
    
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))
    
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080, log_level="info")