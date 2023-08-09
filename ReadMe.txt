You can save your tiff file in "data" and change the "filename" in "demo.m". You can tune the parameters using reference images or set the values for "lam1" and "lam2" directly. The unmixing result will be saved in "result". 

Description: 

Data:
Thirteen reference images of labeled E. coli cells and a mixed biological image

.m files:
demo.m - Demo for unmixing biological fluorescence image data with sparse and low-rank Poisson regression (SL-PRU)
SLPRU.m - Unmixing spectral window with SL-PRU
SLNLS.m - Unmixing spectral window with sparse and low-rank nonnegative least squares
vector_soft_row_w.m - Vectorial soft thresholding operation
Tuning.m - Tuning parameters using reference images
UnmixSpIm.m - Unmixing spectral image with different methods and associated parameter(s)
FigSpIm.m - Displaying estimated abundances for different methods

MIT License

Copyright (c) 2023 Ruogu Wang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
