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
# The seqrepo data is too large to be included in the image
VOLUME /usr/local/share/seqrepo/
VOLUME /home/auto-acmg/seqrepo/master

# ----------------
# AutoACMG runtime
# ----------------

COPY . /auto-acmg

CMD ["pipenv", "run", "uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]

EXPOSE 8000
