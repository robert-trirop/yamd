/*
* Copyright 2021 Robert Schütze
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
* (Created by Robert Schütze on 25.05.2021.)
*/

#include <iostream>
#include "HeaderFiles/verlet.h"
#include "HeaderFiles/types.h"
#include "HeaderFiles/milestones.h"

// this main contains all milestones: uncomment line-wise to get corresponding simulation
int main() {
    MS4(); // uncomment for milestone 4

    // MS5(true, 5, 5, 5); // uncomment for milestone 5 with file output

    // for(int i = 1; i <= 25; i++){MS5(false, 4, 5, i);} // uncomment for milestone 5 time scaling curve

    // MS6(true, 5, 5, 5); // uncomment for milestone 6 with file output

    // for(int i = 1; i <= 50; i++){MS6(false, 4, 5, i);} // uncomment for milestone 6 time scaling curve

    // MS7(); // uncomment for milestone 7

    return 0;
}
