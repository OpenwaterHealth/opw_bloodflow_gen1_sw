'''
Test for control of Thorlabs ITC power supplies
Openwater - Brad Hartl
TO DO's:
Enable logging of laser info to saved scanUI log files
Attempting to overwrite values if failed
How to handle SCPI returns from visaWrite
Need to check if lasers off before turning off TECs?
'''
import pyvisa # Requires running "pip install pyvisa"
import socket
import time

seedLaserThres = [0.000025, 0.00006] # Thresholds (Watts) for value read by "MEAS:POW2?" query [0.00004, 0.00006]
# Note: assumes 1000 mA/W conversion factor on ITC4001
powerUpTimeSeed = 1.5 # Wait period for seed laser to ramp up before power measurement (~1 sec plus OUTP:DEL setting)
powerUpTimeTA = 2 # Wait period for TA
debugLevel = 1 # 1 =  minimal outputs, 2 = all outputs
resourceNameOther = "USB0::0x1313::0x804A::M00543795::0::INSTR" # ITC4002QCL
settingsListSeed = [
    "SOUR2:FUNC TEMP",
    "SOUR2:TEMP 25",
    "OUTP2:STAT 1", # "1" is equivalent to "ON"
    "SOUR:CURR:LIM 0.16",
    "SOUR:CURR 0.15",
    "OUTP:DEL 0.5"]
settingsListTA = [
    "SOUR2:FUNC TEMP",
    "SOUR2:TEMP 25",
    "OUTP2:STAT 1", # "1" is equivalent to "ON"
    "SOUR:FUNC:MODE CURR",
    "SOUR:FUNC:SHAP PULS",
    "TRIG:SOUR EXT",
    "SOUR:PULS:WIDT 0.0001",
    "SOUR:PULS:HOLD WIDT",
    "OUTP:PROT:EXT PROT",
    # "SOUR:PULS:PER 0.067",
    "SOUR:PULS:PER 0.033", # BETATEST
    "SOUR:CURR:LIM 5.1",
    "SOUR:CURR 5",
    # "SOUR:CURR 0.6", # BETATEST
    "OUTP:DEL 1"]

def openLaserConn(): # Opens VISA connection to ITCs
    try:
        resMngr = pyvisa.highlevel.ResourceManager()
        # print(rm.list_resources) # doesn't seem to work
        if socket.gethostname() == "Franklin":
            instrSeed = resMngr.open_resource("USB0::0x1313::0x804A::M00615223::0::INSTR", 0) # VISA for ITC4001 w/o lock
            # instrTA = resMngr.open_resource("USB0::0x1313::0x804A::M00576870::0::INSTR", 0) # VISA for ITC4020 w/o lock
            instrTA = resMngr.open_resource("USB0::0x1313::0x804A::M00613648::0::INSTR", 0) # VISA for ITC4020 w/o lock
            thorLaserStatus = True
            print("Connection to laser drivers opened")
        elif socket.gethostname() == "Schrodinger":
            instrSeed = resMngr.open_resource("USB0::0x1313::0x804A::M00633445::0::INSTR", 0) # VISA for ITC4001 w/o lock
            instrTA = resMngr.open_resource("USB0::0x1313::0x804A::M00582473::0::INSTR", 0) # VISA for ITC4020 w/o lock
            thorLaserStatus = True
            print("Connection to laser drivers opened")
        elif socket.gethostname() == "Strickland":
            instrSeed = resMngr.open_resource("USB0::0x1313::0x804A::M00645524::0::INSTR", 0) # VISA for ITC4001 w/o lock
            instrTA = resMngr.open_resource("USB0::0x1313::0x804A::M00613647::0::INSTR", 0) # VISA for ITC4020 w/o lock
            thorLaserStatus = True
            print("Connection to laser drivers opened")
        else:
            print("No laser drivers available on this scanner")
            instrSeed = instrTA = resMngr = None
            thorLaserStatus = False
    except:
        instrSeed = instrTA = resMngr = None
        thorLaserStatus = False
        print("ERROR: Opening connection to laser drivers failed")
    return instrSeed, instrTA, resMngr, thorLaserStatus

def closeLaserConn(instrSeed, instrTA, resMngr):
    try:
        instrSeed.close()
        instrTA.close()
        resMngr.close()
        print("Connection to laser drivers closed")
    except:
        print("ERROR: Closing connection to laser drivers failed")
        return False
    return True

def turnOnTec(instrSeed, instrTA): # Recall state for each controller, check key values
    try:
        turnOffLaser(instrSeed, instrTA)
        turnOffTec(instrSeed, instrTA)
        outStr1 = instrSeed.query("OUTP:PROT:KEYL:TRIP?") # check key switch
        outStr2 = instrTA.query("OUTP:PROT:KEYL:TRIP?") # check key switch
        if (outStr1[:-1:1] == "1") or (outStr2[:-1:1] == "1"):
            print("ERROR: Keyswitch not enabled")
            return False
        visaWrite(instrSeed, "*RCL 3") # Recall 3 = 4th setting in ITC interface
        visaWrite(instrTA, "*RCL 3") # Recall 3 = 4th setting in ITC interface
        # visaWrite(instrTA, "*RCL 4") # BETATEST
        visaWrite(instrSeed, "OUTP2:STAT ON")
        visaWrite(instrTA, "OUTP2:STAT ON")
        if paramCheck(instrSeed,settingsListSeed,"Seed") and paramCheck(instrTA,settingsListTA,"TA"):
            print("Laser TECs turned on")
        else:
            print("ERROR: Turning on laser TECs failed")
            return False
    except:
        print("ERROR: Turning on laser TECs failed")
        return False
    return True

def turnOffTec(instrSeed, instrTA):
    try:
        if visaWrite(instrSeed, "OUTP2 OFF") and visaWrite(instrTA, "OUTP2 OFF"):
            print("Laser TECs turned off")
        else:
            print("ERROR: Turning off laser TECs failed")
            return False
    except:
        print("ERROR: Turning off laser TECs failed")
        return False
    return True

def turnOnLaser(instrSeed, instrTA):
    try:
        if paramCheck(instrSeed,settingsListSeed,"Seed") and paramCheck(instrTA,settingsListTA,"TA"):
            if visaWrite(instrSeed, "OUTP ON"):
                time.sleep(powerUpTimeSeed)
                outStr = instrSeed.query("MEAS:POW2?")
                if float(outStr[:-1:1]) > seedLaserThres[0] and float(outStr[:-1:1]) < seedLaserThres[1]:
                    print("Success: Seed LD Power is " + outStr[:-1:1])
                    if visaWrite(instrTA, "OUTP ON"):
                        time.sleep(powerUpTimeTA)
                        print("Lasers turned on")
                    else:
                        print("ERROR: Turning on TA failed")
                        return False
                else:
                    print("ERROR: Seed LD Power is " + outStr[:-1:1] + " but expecting between " + str(seedLaserThres[0]) + "and" + str(seedLaserThres[1]))
                    print("ERROR: Turning on lasers failed")
                    return False
            else:
                print("ERROR: Turning on lasers failed")
                return False
        else:
            print("Parmeter checks failed, lasers not turned on")
            return False
    except:
        print("ERROR: Turning on lasers failed")
        return False
    return True

def turnOffLaser(instrSeed, instrTA):
    try:
        if visaWrite(instrTA, "OUTP OFF"):
            time.sleep(0.5) # Waits for half second in case laser is accidentally firing at the time
            if visaWrite(instrSeed, "OUTP OFF"):
                print("Lasers turned off")
            else:
                print("ERROR: Turning off lasers failed")
                return False
        else:
            print("ERROR: Turning off lasers failed")
            return False
    except:
        print("ERROR: Turning off lasers failed")
        return False
    return True

def visaWrite(instr,inputString): # Specifically writing to ITC, immediately checking for errors
    readMessage = instr.write(inputString)
    if debugLevel >= 2:
        if readMessage:
            print("SCPI return: " + str(readMessage))
    errorMessage = instr.query("SYST:ERR?")
    if errorMessage[:-1:1] != "+0,\"No error\"":
        print(errorMessage[:-1:1])
        return False
    return True

def paramCheck(instr,settingsList,instrName): # Checks if crticial settings in list are set correctly
    try:
        for x in settingsList:
            diffTollerance = 0.001 # Tolerance difference acceptable (0.1%)
            splitStr = x.split()
            outStr = instr.query(splitStr[0] + "?")
            if splitStr[1].isalpha():
                paramPass = splitStr[1] == outStr[:-1:1]
            else:
                paramPass = abs(float(splitStr[1]) - float(outStr[:-1:1])) < diffTollerance
            if paramPass:
                if debugLevel >= 2:
                    print("Success: " + splitStr[0] + " set to "  + outStr[:-1:1])
            else:
                print("ERROR: " + splitStr[0] + " set to "  + outStr[:-1:1] + " but expecting " + splitStr[1])
                return False
            time.sleep(0.001)
        print(instrName + " parameter check passed")
    except:
        print("ERROR: " + instrName + " parameter check failed")
        return False
    return True
