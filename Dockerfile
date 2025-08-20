# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.12

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/kg_generator

# Update system and install dependencies.
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    software-properties-common
RUN apt-get install libxrender1 --no-install-recommends -y
RUN apt-get install libxext6 --no-install-recommends -y
# RUN apt-get install libx11-6 --no-install-recommends -y
# RUN apt-get install libglib2.0-0 --no-install-recommends -y
# Copy code and start script (this will place the files in home/username/)
COPY requirements.txt $HOME/kg_generator/requirements.txt
COPY images $HOME/kg_generator/images/
COPY data $HOME/kg_generator/data/
COPY kg_generator.py $HOME/kg_generator/kg_generator.py
# COPY Main.py $HOME/impulse_dashboard/Main.py
COPY kgg_utils.py $HOME/kg_generator/kgg_utils.py
COPY start-script.sh $HOME/kg_generator/start-script.sh

RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8505

HEALTHCHECK CMD curl --fail http://localhost:8505/_stcore/health

ENTRYPOINT ["./start-script.sh"]