# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.11

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/kgg_api

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
# COPY .streamlit $HOME/kgg_api/.streamlit
COPY api_requirements.txt $HOME/kgg_api/api_requirements.txt
# COPY pages $HOME/kgg_api/pages/
COPY images $HOME/kgg_api/images/
COPY data $HOME/kgg_api/data/
COPY api_utils $HOME/kgg_api/api_utils/
# COPY Main.py $HOME/kgg_api/Main.py
COPY kgg_utils.py $HOME/kgg_api/kgg_utils.py
COPY kgg_api.py $HOME/kgg_api/kgg_api.py
COPY start-script.sh $HOME/kgg_api/start-script.sh
COPY .env $HOME/kgg_api/.env

RUN pip install --no-cache-dir -r api_requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8080

ENTRYPOINT ["./start-script.sh"]