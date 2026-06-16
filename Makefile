.PHONY: setup pipeline dashboard clean

setup:
	python -m pip install --upgrade pip
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501 --server.headless true --browser.gatherUsageStats false

clean:
	rm -f cell_counts.db
	rm -rf out