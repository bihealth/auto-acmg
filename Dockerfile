# ==============
# Builc AutoACMG
# ==============

FROM python:3.12-slim
WORKDIR /auto-acmg

# -----------------------
# Install system packages
# -----------------------

RUN apt-get update && apt-get install -y \
    wget \
    gcc \
    libz-dev \
    zlib1g-dev \
    build-essential \
    tabix \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

RUN pip install pipenv

# ----------------------
# Set up the environment
# ----------------------

COPY Pipfile Pipfile.lock .

RUN PIPENV_VENV_IN_PROJECT=1 pipenv install --deploy

# Fix the seqrepo bug https://github.com/biocommons/biocommons.seqrepo/issues/166
RUN sed -i -e 's/if aliases_cur.fetchone() is not None/if next(aliases_cur, None) is not None/' \
    .venv/lib/python3.12/site-packages/biocommons/seqrepo/cli.py

# Set up seqrepo
# RUN mkdir -p /usr/local/share/seqrepo \
#     && chown -R root:root /usr/local/share/seqrepo
# RUN pipenv run seqrepo init -i auto-acmg
# RUN pipenv run seqrepo fetch-load -i auto-acmg -n RefSeq NC_000001.10 NC_000002.11 NC_000003.11 \
#     NC_000004.11 NC_000005.9 NC_000006.11 NC_000007.13 NC_000008.10 NC_000009.11 NC_000010.10 \
#     NC_000011.9 NC_000012.11 NC_000013.10 NC_000014.8 NC_000015.9 NC_000016.9 NC_000017.10 \
#     NC_000018.9 NC_000019.9 NC_000020.10 NC_000021.8 NC_000022.10 NC_000023.10 NC_000024.9 \
#     NC_012920.1 NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 \
#     NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 \
#     NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 \
#     NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10 NC_012920.1

# The seqrepo data is too large to be included in the image
VOLUME /usr/local/share/seqrepo/
VOLUME /home/auto-acmg/seqrepo/master

# ----------------
# AutoACMG runtime
# ----------------

COPY . /auto-acmg

CMD ["pipenv", "run", "uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]

EXPOSE 8000
