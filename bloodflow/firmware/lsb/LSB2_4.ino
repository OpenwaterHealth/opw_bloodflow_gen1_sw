//Laser Safety Board
// By Peter Bradbury
// Rev 2

#define DIP 5
#define COIL 10
#define LED A5
#define PWR A2
#define Bled A3
#define HRTBT_OUT 11
#define HRTBT_IN A4
#define SR_FBK A0



int PD = 2;
int PDs = 1;
int DIPs = 0;
int LEDstate = 0;
int MCINs = 1;
int flag = 0;

unsigned long duration;
unsigned long duration2;
unsigned long previousMillis = 0;
const long interval = 1000;
int PWRraw;
float PWRv;
unsigned long startMicros;
unsigned long endMicros;
unsigned long totalMicros;
//unsigned long startMicros2;
//unsigned long totalMicros2;
float SR_FBKv;

int flag2 = 0;
int flag3 = 1;
unsigned long HRTBT_duration;

unsigned long timer = 0;
unsigned long endTimer = 0;

void setup() {
  pinMode(DIP, INPUT);
  pinMode(COIL, OUTPUT);
  pinMode(LED, OUTPUT);
  pinMode(PD, INPUT);
  pinMode(Bled, OUTPUT);
  //pinMode(PWR, INPUT);
    pinMode(HRTBT_OUT, OUTPUT);
  pinMode(HRTBT_IN, INPUT);
    pinMode(SR_FBK, INPUT);

  digitalWrite(COIL, HIGH);
  Serial1.begin(115200);
  //delay(1000);
//    while (!Serial)
//  {
//    ; // wait for serial port to connect. Needed for Leonardo only
//  }
  Serial1.println("setup");
analogWrite(HRTBT_OUT, 125);
}

void loop() {

//DIPs = digitalRead(DIP);
//if (DIPs == LOW){
// analogWrite(MCOUT, 125); 
//}
//
//if (DIPs == HIGH){
//  analogWrite(MCOUT, 255);
//}

//  duration = pulseIn(PD, HIGH);
//  Serial.print("Pulse Width: ");
//  Serial.print(duration);
//  Serial.println(" uS");

//startMicros2 = micros();
//  duration2 = pulseIn(MCIN, HIGH, 800);
//  Serial.print("Pulse Width: ");
//  Serial.print(duration2);
//  Serial.println(" uS");

//  startMicros = micros();
//MCINs = digitalRead(MCIN);
PDs = digitalRead(PD);
if (flag == 0){
if (PDs == LOW){
  startMicros = micros();
  flag = 1;
}
//  endMicros = micros();
//  //totalMicros = endMicros - startMicros;
//  totalMicros2 = endMicros - startMicros2;
//  //Serial.print("total micros:");
//  //Serial.println(totalMicros);
//  Serial.print("total micros 2:  ");
//  Serial.println(totalMicros2);
//  Serial.println("LASER KILL");
//  delay(1000);
}
PDs = digitalRead(PD);
if(PDs == HIGH){
  startMicros = micros();
}
if (flag == 1){
  endMicros = micros();
  totalMicros = endMicros - startMicros;
  if (totalMicros >= 500){
    digitalWrite (COIL, LOW);
    Serial1.println("LOW");
    flag = 0;
    delay(5000);
  }
}
//  digitalWrite(Bled, HIGH);
//  delay(1000);
//  digitalWrite(Bled, LOW);

  
//  PWRraw = analogRead(PWR);
//  delay(50);
//  PWRraw = analogRead(PWR);
//  delay(50);
//  PWRv = PWRraw*(5.0/1023.0);
//   Serial.print("PWR voltage: ");
//  Serial.print(PWRv);
//  Serial.println(" V");
  
  unsigned long currentMillis = millis();

 if (currentMillis - previousMillis >= interval) {
    // save the last time you blinked the LED
    previousMillis = currentMillis;

     //if the LED is off turn it on and vice-versa:
    if (LEDstate == LOW) {
      LEDstate = HIGH;
    } else {
      LEDstate = LOW;
    }

//     set the LED with the ledState of the variable:
    digitalWrite(LED, LEDstate);
    digitalWrite(Bled, LEDstate);
      }
if (duration == 0){
  if (flag2 == 0){
  timer = millis();
  flag2 = 1;
  }
  endTimer = millis();
  if (endTimer - timer >= 200000){
    flag2 = 0;
    flag3 = 1;
  }
}

if (flag3 == 1){
  HRTBT_duration = pulseIn(HRTBT_IN, HIGH);
  SR_FBKv = analogRead(SR_FBK);
  flag3 = 0;
}

if (HRTBT_duration == 0){
  digitalWrite (COIL, LOW);
}

if (COIL == LOW && SR_FBKv >= 3.5){
  Serial.println("MCU SR control malfunciton or SR malfunction");
  Serial.println("MCU requesting SR coil low, but SR coil is high or contacts welded");
}

if (COIL == HIGH && SR_FBKv <= 1.5){
    Serial.println("MCU SR control malfunciton or SR malfunction");
  Serial.println("MCU requesting SR coil HIGH, but SR coil is LOW or contacts welded");
}


  
// DIPs = digitalRead(DIP);
//if (DIPs == HIGH){
//  digitalWrite (COIL, HIGH);
//  Serial.println("high");
//  
//}
//if (duration >= 490 && duration <= 510){
//  digitalWrite (COIL, HIGH);
//    Serial.println("HIGH");
//}
//if (duration <= 489){ 
//  digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
//if (duration >= 512){
//   digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
}
