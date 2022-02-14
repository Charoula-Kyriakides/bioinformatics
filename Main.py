#from functions import *

# def fibonachi(n):
#     if n==0:
#         return 0
#     if n==1:
#         return 1
#     a = 0
#     b = 1
#     fib = 0
#     for _ in range(n-1):
#         fib = a+b
#         a = b
#         b = fib
#     return fib


import pytesseract
pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR\tesseract'
print(pytesseract.image_to_string(r'D:\examplepdf2image.png'))