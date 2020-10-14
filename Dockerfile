# Set the base image using miniconda 
FROM continuumio/miniconda3:4.3.27

# Set environmental variable(s)
ENV ACCEPT_INTEL_PYTHON_EULA=yes

# Set working directory 
WORKDIR /home/notebooks

# Add requirements file 
ADD requirements.txt /app/

# Installs, clean, and update   
RUN pip install -r /app/requirements.txt

CMD jupyter notebook --port=8888 --ip=0.0.0.0 --allow-root

    