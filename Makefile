.PHONY: test lint demo build validate serve clean

test:
	python -m pytest -q

lint:
	ruff check pipeline tests
	npx --yes prettier@3.9.6 --check index.html assets config

demo:
	python -m pipeline.build --output data/v3 --max-records 750

build:
	python -m pipeline.build --output data/v3 --allow-partial

validate:
	python -m pipeline.validate data/v3

serve:
	python -m http.server 8000

clean:
	rm -rf .cache build _site
