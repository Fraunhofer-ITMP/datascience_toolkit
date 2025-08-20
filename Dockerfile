# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.12

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/drug_likeness_assessment

# Update system and install dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    chromium \
    libnss3 \
    libasound2 \
    fonts-liberation \
    libatk-bridge2.0-0 \
    libgtk-3-0 \
    libx11-xcb1 \
    && rm -rf /var/lib/apt/lists/*
# RUN apt-get install libx11-6 --no-install-recommends -y
# RUN apt-get install libglib2.0-0 --no-install-recommends -y
# Copy code and start script (this will place the files in home/username/)
# COPY .streamlit $HOME/drug_likeness_assessment/.streamlit
COPY requirements.txt $HOME/drug_likeness_assessment/requirements.txt
# COPY pages $HOME/drug_likeness_assessment/pages/
COPY images $HOME/drug_likeness_assessment/images/
COPY data $HOME/drug_likeness_assessment/data/
COPY drug_likeness_assessment.py $HOME/drug_likeness_assessment/drug_likeness_assessment.py
# COPY Main.py $HOME/drug_likeness_assessment/Main.py
# COPY kgg_utils.py $HOME/drug_likeness_assessment/kgg_utils.py
COPY start-script.sh $HOME/drug_likeness_assessment/start-script.sh

RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8503

HEALTHCHECK CMD curl --fail http://localhost:8503/_stcore/health

ENTRYPOINT ["./start-script.sh"]