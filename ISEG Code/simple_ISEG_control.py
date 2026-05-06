# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 14:39:37 2026

@author: iontrap
"""

from mpod_control import MPOD
import time
import numpy as np

mp = MPOD()

slot = 0

#l1
channel1 = 0
#L2
channel2 = 1
#L3
channel3 = 2
#L4
channel4 = 3
#L5
channel5 = 4
#R1
channel6 = 5
#R2
channel7 = 6
#R3
channel8 = 7
#R4
channel9 = 8
#R5
channel10 = 9

L1 = 1.0
L2 = 0.0
L3 = 0.0
L4 = 0.0
L5 = 0.0
R1 = 0.0
R2 = 0.0
R3 = 0.0
R4 = 0.0
R5 = 0.0



mp.set_voltage(slot, channel1, L1)

mp.set_voltage(slot, channel2, L2)

mp.set_voltage(slot, channel3, L3)

mp.set_voltage(slot, channel4, L4)

mp.set_voltage(slot, channel5, L5)

mp.set_voltage(slot, channel6, R1)

mp.set_voltage(slot, channel7, R2)

mp.set_voltage(slot, channel8, R3)

mp.set_voltage(slot, channel9, R4)

mp.set_voltage(slot, channel10, R5)

mp.on(slot, channel1)
mp.on(slot, channel2)
mp.on(slot, channel3)
mp.on(slot, channel4)
mp.on(slot, channel5)
mp.on(slot, channel6)
mp.on(slot, channel7)
mp.on(slot, channel8)
mp.on(slot, channel9)
mp.on(slot, channel10)

input("\nPress ENTER to turn OFF all channels...")

# Turn OFF all channels
for ch in range(10):
    mp.off(slot, ch)

print("All channels OFF.")