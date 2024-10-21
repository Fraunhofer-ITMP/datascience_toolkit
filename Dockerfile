# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.10-slim

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/datascience_toolkit

# Update system and install dependencies.
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    software-properties-common
RUN apt-get install libxrender1 --no-install-recommends -y
RUN apt-get install libxext6 --no-install-recommends -y
# RUN apt-get install libx11-6 --no-install-recommends -y
# RUN apt-get install libglib2.0-0 --no-install-recommends -y

# Copy code and start script (this will place the files in home/username/)
COPY .streamlit $HOME/datascience_toolkit/.streamlit
COPY requirements.txt $HOME/datascience_toolkit/requirements.txt
COPY pages $HOME/datascience_toolkit/pages/
COPY images $HOME/datascience_toolkit/images/
COPY data $HOME/datascience_toolkit/data/
COPY Main.py $HOME/datascience_toolkit/Main.py
COPY kgg_utils.py $HOME/datascience_toolkit/kgg_utils.py
COPY start-script.sh $HOME/datascience_toolkit/start-script.sh

RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["./start-script.sh"]