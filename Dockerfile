# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.12

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/impulse_dashboard

# Update system and install dependencies.
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    software-properties-common
RUN apt-get install libxrender1 --no-install-recommends -y
RUN apt-get install libxext6 --no-install-recommends -y
# RUN apt-get install libx11-6 --no-install-recommends -y
# RUN apt-get install libglib2.0-0 --no-install-recommends -y
# Copy code and start script (this will place the files in home/username/)
# COPY .streamlit $HOME/impulse_dashboard/.streamlit
COPY requirements.txt $HOME/impulse_dashboard/requirements.txt
# COPY pages $HOME/impulse_dashboard/pages/
COPY images $HOME/impulse_dashboard/images/
COPY data $HOME/impulse_dashboard/data/
COPY 6_OpenScreen_Impulse.py $HOME/impulse_dashboard/6_OpenScreen_Impulse.py
# COPY Main.py $HOME/impulse_dashboard/Main.py
# COPY kgg_utils.py $HOME/impulse_dashboard/kgg_utils.py
COPY start-script.sh $HOME/impulse_dashboard/start-script.sh

RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["./start-script.sh"]