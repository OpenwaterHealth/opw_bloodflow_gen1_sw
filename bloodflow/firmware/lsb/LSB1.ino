//Laser Safety Board
// By Peter Bradbury
// Rev 1

#define DIP 5
#define COIL 10
#define LED A5
#define PWR A2
#define BLED A3


int PD = 2;
int DIPs = 0;
int LEDstate = 0;
int flag1 = 0;
unsigned long duration;
unsigned long previousMillis = 0;
const long interval = 1000;
int PWRraw;
float PWRv;

void setup() {
  pinMode(DIP, INPUT);
  pinMode(COIL, OUTPUT);
  pinMode(LED, OUTPUT);
  pinMode(PD, INPUT);
  pinMode(BLED, OUTPUT);
  //pinMode(PWR, INPUT);
  
  Serial1.begin(115200);
  delay(10);
  digitalWrite(COIL, HIGH);
//delay(10000);
//    while (!Serial)
//  {
//    ; // wait for serial port to connect. Needed for Leonardo only
//  }
  Serial.println("setup");

}

void loop() {

  duration = pulseIn(PD, LOW);
  Serial.print("Pulse Width: ");
  Serial.print(duration);
  Serial.println(" uS");
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

    // if the LED is off turn it on and vice-versa:
    if (LEDstate == LOW) {
      LEDstate = HIGH;
    } else {
      LEDstate = LOW;
    }

    // set the LED with the ledState of the variable:
    digitalWrite(LED, LEDstate);
  
  }
  
 //DIPs = digitalRead(DIP);
//if (DIPs == HIGH){
//  digitalWrite (COIL, HIGH);
//  Serial.println("high");
//  
//}
//
//if (DIPs == LOW){
//  digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
if (duration <= 430){
  digitalWrite (COIL, HIGH);
    //Serial.println("HIGH");
    Serial.print ("HIGH");
}


//if (duration <= 310 && duration >= 150){ 
//  digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
if (duration >= 440){
   digitalWrite (COIL, LOW);
  Serial.println("LOW");
}
}
