pyinstaller --noconfirm --onefile --windowed --distpath ./ --icon ./Resources/shadow.ico --clean Shadow.py

rm Shadow.spec
rm -r __pycache__
rm -r ./GUI/__pycache__
rm -r build
