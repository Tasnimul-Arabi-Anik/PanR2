FROM mambaorg/micromamba:1.5.10

ARG SETUP_ABRICATE_DB=false

WORKDIR /opt/PanR2

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

COPY . /opt/PanR2
RUN /opt/conda/bin/pip install -e .

RUN if [ "$SETUP_ABRICATE_DB" = "true" ]; then /opt/conda/bin/abricate --setupdb; fi

ENV PATH="/opt/conda/bin:${PATH}"
ENTRYPOINT ["panr"]
CMD ["--help"]
