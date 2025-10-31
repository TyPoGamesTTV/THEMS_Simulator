@echo off
python -m venv venv
call venv\Scripts\activate
pip install -r requirements.txt
python thems_sim.py --beta 110
echo Done! Check thems_simulation.png
pause