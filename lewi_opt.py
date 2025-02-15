# Required imports

import datetime
import math
import matplotlib.pyplot
import mealpy
import numpy
import os
import random
import scipy.optimize
import subprocess
import threading
import time
import sys

# Global constants

# coil.geo

numberOfTurns = 3
numberOfTurns2 = 8
numberOfTurns3 = 3
numberOfTurns4 = 8
numberOfUpperTurns = 3
# radius = 0.03
turnDist = math.pi / 10

# levitation_melting_2d.geo
surfaceTension = 1.557
charge_R = 0.1
chargeY = 0
conductorRadius = 0.003
margin = 0.003
turn_r = 0.003

# Z daszkiem
# -2.7 dla 18 mm 10 kHz
# -1.45 dla 15 mm 10 kHz
# -0.5 dla 12 mm 10 kHZ

charge_r = 0.012 # ROZMIAR WSADU
forceSurfaceTension = -1.45 # TU STOSOWAĆ Z POWYŻSZEJ TABELKI DLA POSZCZEGÓLNYCH ROZMIARÓW WSADU

# other
csvPath = 'opt.csv'
logPath = 'opt.log'
logColumnDelimiter = ', '
logRowDelimiter = '\n'

# Threading
GLOBAL_LOCK = threading.Lock()

# Optimization progress chart
optimalY_totalForce = None
optimalArguments = None

# def distance(argument):
#   coilX = argument[0]
#   coilY = argument[1]
#   turnStart = argument[2]
#   distance = math.sqrt((coilX - chargeX) ** 2 + (coilY - chargeY) ** 2)
#   return distance

lastCriterionValue = 0
requiredCriterionTolerance = 0.001
currentBestCriterion = None
totalDuration = 0

def reset():
  global currentBestCriterion, totalDuration
  currentBestCriterion = None
  totalDuration = 0

def log(message):
  with open(logPath, 'a') as logFile:
    logFile.write(f'[{datetime.datetime.now().isoformat()}Z] ')
    logFile.write(message)
    logFile.write(logRowDelimiter)
    logFile.flush()
    print(message)

def csv(*args):
  with open(csvPath, 'a') as csvFile:
    csvFile.write(logColumnDelimiter.join(str(value) for value in args))
    csvFile.write(logRowDelimiter)
    csvFile.flush()

def readTotalForce():
  # Output (GlobalForce.txt)
  with open('GlobalForce_r.txt', 'r') as forceFile:
    force = forceFile.read().split()
    #forceX = float(force[1])

  with open('GlobalForce.txt', 'r') as forceFile:
    force = forceFile.read().split()
    forceX = float(force[1])
    forceY = float(force[2])

  # Output (GlobalGravity.txt)
  with open('GlobalGravity.txt', 'r') as gravityFile:
    gravity = gravityFile.read().split()
    gravity = float(gravity[1])

  forceComponentX = (forceX - forceSurfaceTension) ** 2
  forceComponentY = (forceY - gravity) ** 2
  totalForce = math.sqrt(forceComponentX + forceComponentY)
  return totalForce

# Argument constraints
# C1. d < R - (r_ch + r_t + margin)

def turnPositions():
  number_of_turns = 8
  turn_start_cos = -math.pi / 2
  turn_start_sin = -math.pi / 2
  coil_x = 0.1
  coil_y = 0.03
  radius = 0.06
  turn_dist=math.pi / 20
  indx = []
  indy = []

  for i in range(number_of_turns):
    indx.append( coil_x + radius * math.cos(turn_start_cos + turn_dist * (i - (number_of_turns - 1) / 2)))
    indy.append( coil_y + radius * math.sin(turn_start_sin + turn_dist * (i - (number_of_turns - 1) / 2)))

  xs = indx[0]
  ys = indy[0]
  xe = indx[len(indx) - 1]
  ye = indy[len(indy) - 1]

  for i in range(numberOfUpperTurns):
    indx.insert(0, xs)
    indy.insert(0, ys + 0.010 * (i + 1))
    indx.append(xe)
    indy.append(ye + 0.010 * (i + 1))

  ys = max(indy)

  for i in range(numberOfTurns4):
    indy.append(ys + 0.02)
    indx.append(xs + (xe - xs) / (numberOfTurns4 - 1) * i)

  # for i in range(numberOfUpperTurns):
  #   indx.insert(0, xs)
  #   indy.insert(0, ys + 0.010 * (i + 6))
  #   indx.append(xe)
  #   indy.append(ye + 0.010 * (i + 6))

  return (indx, indy)

def funC1(argument):
  global charge_r, conductorRadius, margin
  coilX = argument[0]
  coilY = argument[1]
  turnStartCos = argument[2]
  turnStartSin = argument[3]
  radius = argument[4]
  # coilY = argument[1] + radius

  if numberOfTurns2 > 0:
    coilX2 = argument[5]
    coilY2 = argument[6]
    turnStartCos2 = argument[7]
    turnStartSin2 = argument[8]
    radius2 = argument[9]
    # coilY2 = argument[6] + radius2

  (indx, indy) = turnPositions(numberOfTurns, turnStartCos, turnStartSin, coilX, coilY, radius, turnDist * 0.03 / radius)
  minvalue = 0

  for i in range(len(indx)):
    turnChargeDistance = math.sqrt((charge_R - indx[i]) ** 2 + (chargeY - indy[i]) ** 2) - (charge_r + conductorRadius + margin)

    if turnChargeDistance < 0:
      minvalue += turnChargeDistance

    turnAxisDistance = indx[i] - conductorRadius - margin

    if turnAxisDistance < 0:
      minvalue += turnAxisDistance

    for j in range (i):
      turnTurnDistance = math.sqrt((indx[j] - indx[i]) ** 2 + (indy[j] - indy[i]) ** 2) - (2 * conductorRadius + 0.001)

      if turnTurnDistance < 0:
        minvalue += turnTurnDistance

  coilLeftDistance = (charge_R - charge_r - turn_r) - min(indx)

  if coilLeftDistance < 0:
    minvalue += coilLeftDistance

  coilRightDistance = max(indx) - (charge_R + charge_r)

  if coilRightDistance < 0:
    minvalue += coilRightDistance

  coil1chargeY = chargeY - min(indy)

  if coil1chargeY < 0:
    minvalue += coil1chargeY

  if numberOfTurns2>0:
    (indx2, indy2) = turnPositions(numberOfTurns2, turnStartCos2, turnStartSin2, coilX2, coilY2, radius2, turnDist * 0.03 / radius2)

    for i in range(numberOfTurns2):
      turnChargeDistance = math.sqrt((charge_R - indx2[i]) ** 2 + (chargeY - indy2[i]) ** 2) - (charge_r + conductorRadius + margin)

      if turnChargeDistance < 0:
        minvalue += turnChargeDistance

      turnAxisDistance = indx2[i] - conductorRadius - margin

      if turnAxisDistance < 0:
        minvalue += turnAxisDistance

      for j in range (i):
        turnTurnDistance = math.sqrt((indx2[j] - indx2[i]) ** 2 + (indy2[j] - indy2[i]) ** 2) - (2 * conductorRadius + 0.001)

        if turnTurnDistance < 0:
          minvalue += turnTurnDistance

    coilLeftDistance = min(indx2) - min(indx)

    if coilLeftDistance < 0:
      minvalue += coilLeftDistance

    coilRightDistance = max(indx) - max(indx2)

    if coilRightDistance < 0:
      minvalue += coilRightDistance

    coil2chargeY = max(indy2) - chargeY

    if coil2chargeY < 0:
      minvalue += coil2chargeY

    for i1 in range(len(indx)):
      for i2 in range(numberOfTurns2):
        turnTurnDistance = math.sqrt((indx2[i2] - indx[i1]) ** 2 + (indy2[i2] - indy[i1]) ** 2) - (2 * conductorRadius + 0.001)

        if turnTurnDistance < 0:
          minvalue += turnTurnDistance

  # distance = math.sqrt((coilX - charge_R) ** 2 + (coilY - chargeY) ** 2)
  # value1 = (radius - (charge_r + conductorRadius + margin) - distance)
  # value2 = (coilX - radius - conductorRadius)
  # if value1 < value2:
  #   return value1
  # else:
  #   return value2

  return minvalue

constraintC1 = {
  'type': 'ineq',
  'fun': funC1
}

argumentConstraints = [constraintC1]

class StopOptimization(Exception):
  pass

def outputCallback(result, duration):
  print(result, duration)

def runGMSHandGetDP(argument):
  with GLOBAL_LOCK:
    # Measure iteration time
    timeStart = time.time()

    # Optimized variables
    global optimalArguments, optimalY_totalForce, currentBestCriterion, totalDuration

    current1 = argument[0]
    current2 = argument[1]
    current3 = argument[2]
    current4 = argument[3]
    (indx, indy) = turnPositions()

    # Input (coil.geo)
    with open('coil.geo', 'w') as coilFile:
      # coilFile.write(f'coil_x = {coilX};\n')
      # coilFile.write(f'coil_y = {coilY};\n')
      coilFile.write(f'number_of_turns = {numberOfTurns};\n')
      coilFile.write(f'number_of_turns2 = {numberOfTurns2};\n')
      coilFile.write(f'number_of_turns3 = {numberOfTurns3};\n')
      coilFile.write(f'number_of_turns4 = {numberOfTurns4};\n')
      coilFile.write(f"indx() = {str(indx).replace('[', '{').replace(']', '}')};\n")
      coilFile.write(f"indy() = {str(indy).replace('[', '{').replace(']', '}')};\n")
      coilFile.write(f'charge_r = {charge_r};\n')
      coilFile.write(f'charge_R = {charge_R};\n')
      coilFile.write(f"Current1 = {current1};\n")
      coilFile.write(f"Current2 = {current2};\n")
      coilFile.write(f"Current3 = {current3};\n")
      coilFile.write(f"Current4 = {-current4};\n")

      # turnDistScale = turnDist * 0.03 / radius
      # coilFile.write(f'turn_dist = {turnDistScale};\n')
      # coilFile.write(f'turn_start = {turnStart};\n')
      coilFile.flush()

    # Run GMSH and wait
    # 31.07.2024 - Using 'subprocess.run' instead of 'os.system' for calling external program
    # os.system('./gmsh -v 0 -3 levitation_melting_2d.geo')

    while True:
      ok = 1

      try:
        subprocess.run(['./gmsh', '-v', '0', '-3', 'levitation_melting_2d.geo'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, timeout = 30)

      except subprocess.TimeoutExpired:
        ok = 0

      else:
        try:
          subprocess.run(['./getdp', '-v', '0', 'levitation_melting_2d.pro', '-solve', 'Melting', '-pos', 'full', '-v2'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, timeout = 30)

        except subprocess.TimeoutExpired:
          ok = 0

        else:
          if ok == 1:
            break

    totalForce = readTotalForce()

    if optimalY_totalForce == None or optimalY_totalForce > totalForce:
      optimalY_totalForce = totalForce
      optimalArguments = argument

    timeEnd = time.time()
    duration = timeEnd - timeStart
    totalDuration = totalDuration + duration

    if currentBestCriterion == None or totalForce < currentBestCriterion:
      # Znaleziono najlepszą historycznie wartość
      outputCallback(result = totalForce, duration = totalDuration)

    currentBestCriterion = totalForce if currentBestCriterion == None or totalForce < currentBestCriterion else currentBestCriterion
    log(f'Run runGMSHandGetDP() c1 = {current1}, c2 = {current2}, c3 = {current3}, c4 = {current4}, result = {totalForce}, currentBestCriterion = {currentBestCriterion}, duration = {duration} [s], totalDuration = {totalDuration} [s]')
    csv(f'{datetime.datetime.now().isoformat()}Z', current1, current2, current3, current4, totalForce, currentBestCriterion, duration, totalDuration)

    if totalForce < requiredCriterionTolerance:
      # with open('D:\\SGolak\\Lewitacja\\lewitacja\\lastOptimal.txt', 'w') as file:
      #   file.write(' '.join(map(str, optimalArguments)))

      raise StopOptimization()

    return totalForce

initialCurrent1 = 100
initialCurrent2 = 100
initialCurrent3 = 100
initialCurrent4 = 100
initialArgument = [initialCurrent1, initialCurrent2, initialCurrent3, initialCurrent4]
# initialArgument = [414.95834185, 576.983245, 215.77349556];
# with open('D:\\SGolak\\Lewitacja\\lewitacja\\lastOptimal.txt', 'r') as file:
#   content = file.read()
#   initialArgument = list(map(float, content.split()))

boundsCurrent1 = (100, 2000)
boundsCurrent2 = (100, 2000)
boundsCurrent3 = (100, 2000)
boundsCurrent4 = (100, 2000)
# boundsCurrent1 = (220, 230)
# boundsCurrent2 = (630, 640)
# boundsCurrent3 = (60, 70)
argumentBounds = (boundsCurrent1, boundsCurrent2, boundsCurrent3, boundsCurrent4)

def stop_optimization(xk, convergence = None):
  global lastCriterionValue

  if lastCriterionValue < requiredCriterionTolerance:
    raise StopOptimization(f"Value of the objective function below {requiredCriterionTolerance}")

# Optimization

if os.path.isfile(logPath) == True:
  os.remove(logPath)

if os.path.isfile(csvPath) == True:
  os.remove(csvPath)

def callback(result, duration):
  global X, Y
  X = numpy.append(X, duration / 3600)
  Y = numpy.append(Y, result)
  print(f'FOUND BEST ({result}) after {duration} [s]')

# Callback function for logging
outputCallback = callback

for CHARGE_R, FORCE_SURFACE_TENSION in [
  [0.012, -0.5],
  [0.015, -1.45],
  [0.018, -2.7]
]:
  charge_r = CHARGE_R
  forceSurfaceTension = FORCE_SURFACE_TENSION
  log(f'*** Optymalizacja wariantu parametrów dla charge_r = {charge_r}, forceSurfaceTension = {forceSurfaceTension} ***')
  csv(charge_r, forceSurfaceTension)
  matplotlib.pyplot.figure()

  # Optymalizacja dla algorytmów NELDER-MEAD, POWELL, COBYLA
  for ALGORITHM in ['Nelder-Mead', 'Powell', 'COBYLA']:
    log(f'*** Optymalizacja wariantu dla algorytmu "{ALGORITHM}" ***')
    csv(ALGORITHM)
    reset()

    try:
      X = numpy.array([])
      Y = numpy.array([])
      result = scipy.optimize.minimize(runGMSHandGetDP, initialArgument, bounds = argumentBounds, method = ALGORITHM)

    except StopOptimization:
      log('*** Optymalizacja odnalazla zadawalajace minimum ***')
      matplotlib.pyplot.plot(X, Y, label = ALGORITHM)

  # Optymalizacja dla algorytmu PSO
  try:
    log(f'*** Optymalizacja wariantu dla algorytmu "PSO (Particle Swarm Optimization)" ***')
    csv('PSO')
    reset()
    X = numpy.array([])
    Y = numpy.array([])
    optimizer = mealpy.swarm_based.PSO.OriginalPSO()
    optimizer.solve({'obj_func': runGMSHandGetDP, 'bounds': mealpy.FloatVar(lb = numpy.ones(4) * 100, ub = numpy.ones(4) * 2000), 'minmax': 'min'})

  except StopOptimization:
    log('*** Optymalizacja odnalazla zadawalajace minimum ***')
    matplotlib.pyplot.plot(X, Y, label = 'PSO')

  # Optymalizacja dla algorytmu GWO
  try:
    log(f'*** Optymalizacja wariantu dla algorytmu "GWO (Grey Wolf Optimization)" ***')
    csv('GWO')
    reset()
    X = numpy.array([])
    Y = numpy.array([])
    optimizer = mealpy.swarm_based.GWO.OriginalGWO()
    optimizer.solve({'obj_func': runGMSHandGetDP, 'bounds': mealpy.FloatVar(lb = numpy.ones(4) * 100, ub = numpy.ones(4) * 2000), 'minmax': 'min'})

  except StopOptimization:
    log('*** Optymalizacja odnalazla zadawalajace minimum ***')
    matplotlib.pyplot.plot(X, Y, label = 'GWO')

  # Optymalizacja dla algorytmu GA
  try:
    log(f'*** Optymalizacja wariantu dla algorytmu "GA (Genetic Algorithm)" ***')
    csv('GA')
    reset()
    X = numpy.array([])
    Y = numpy.array([])
    optimizer = mealpy.evolutionary_based.GA.BaseGA()
    optimizer.solve({'obj_func': runGMSHandGetDP, 'bounds': mealpy.FloatVar(lb = numpy.ones(4) * 100, ub = numpy.ones(4) * 2000), 'minmax': 'min'})

  except StopOptimization:
    log('*** Optymalizacja odnalazla zadawalajace minimum ***')
    matplotlib.pyplot.plot(X, Y, label = 'GA')

  matplotlib.pyplot.grid()
  matplotlib.pyplot.legend(loc = 'upper right')
  matplotlib.pyplot.title(f'{CHARGE_R * 1000} [mm]')
  matplotlib.pyplot.xlabel('Execution time [h]')
  matplotlib.pyplot.xscale('log')
  matplotlib.pyplot.ylabel('Total criterion [?]')
  matplotlib.pyplot.savefig(f'{CHARGE_R * 1000}mm.png')

# print(f'Optimal current1: {optimalArguments[0]}')
# print(f'Optimal current2: {optimalArguments[1]}')
# print(f'Optimal current3: {optimalArguments[2]}')
# print(f'Optimal current4: {optimalArguments[3]}')
# print(f'Optimal totalForce: {optimalY_totalForce}')
# print(runGMSHandGetDP(optimalArguments))

# with open('D:\\SGolak\\Lewitacja\\lewitacja\\lastOptimal.txt', 'w') as file:
#   file.write(' '.join(map(str, optimalArguments)))
