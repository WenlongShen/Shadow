pyinstaller --noconfirm --onefile --windowed --distpath ./ --icon ./resources/shadow.ico --clean shadow.py

rm shadow.spec
rm -r __pycache__
rm -r build
