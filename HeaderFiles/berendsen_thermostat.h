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
* (Created by Robert Schütze on 17.06.2021.)
*/


#ifndef YAMD_BERENDSEN_THERMOSTAT_H
#define YAMD_BERENDSEN_THERMOSTAT_H

#include "atoms.h"
#include "constants.h"
#include "useful_functions.h"

// inputs: target temperature T0, timestep dt and relaxation time tau
void berendsen_thermostat(Atoms &atoms, double T0, double dt, double tau);

#endif //YAMD_BERENDSEN_THERMOSTAT_H
