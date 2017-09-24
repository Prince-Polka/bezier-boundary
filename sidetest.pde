/*
 approximate cuberoot, 
 ported from c code @
 http://www.hackersdelight.org/hdcodetxt/acbrt.c.txt
 */
final static float acbrt(float x) {
  int xi = Float.floatToIntBits(x);
  xi = (xi>>2)+(xi>>4);
  xi += xi>>4;
  xi += xi>>8;
  xi += 709965728;
  float xf = Float.intBitsToFloat(xi);
  xf = 0.33333333f*(xf+xf+x/(xf*xf));  // Newton step.
  return xf;
}
/* 
 approximate arctan , 
 removed if (y < 0) as y is always positive
 original code from njuffa @ https://math.stackexchange.com/a/1105038 
 */
final static float fastan2(float y, float x) {
  float aX = Math.abs(x);
  float a, s, r;
  a=Math.min(aX, y)/Math.max(aX, y); 
  s = a*a;
  r = ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
  if (y>aX)r=HALF_PI-r;
  if (x<0)r=PI-r;
  return r;
}
/*
approximate modulo
 suprinsignly faster than the & operator for floating point types
 */

final static float fastmod(float a, float b) {
  return a - b*(long)(a/b);
}
final static float 
  invfact3 = - 1 / 6f, 
  invfact5 =   1 / 120f, 
  invfact7 = - 1 / 5040f, 
  invfact9 =   1 / 362880f, 
  invfact11 = - 1 / 39916800f, 
  invfact13 =   1 / 6227020800f, 
  invfact15 = - 1 / 1307674368000f, 
  invfact17 =   1 / 355687428096000f, 
  invfact19 = - 1 / 121645100408832000f, 

  PI = (float)Math.PI, 
  TWO_PI =(float)Math.PI*2, 
  HALF_PI = (float)Math.PI/2, 
  SIXTH_PI = (float)Math.PI/6, 
  CBRT2 = 1.2599210498948731f;
;
final static float fastsin(float x) {
  final float x0 = fastmod(x, TWO_PI); 
  x = fastmod(x, PI); 
  final float x2=x*x;
  x = x 
    + (x*=x2) * invfact3
    + (x*=x2) * invfact5
    + (x*=x2) * invfact7
    + (x*=x2) * invfact9
    + (x*=x2) * invfact11
    + (x*=x2) * invfact13
    + (x*=x2) * invfact15
    + (x*=x2) * invfact17
    + (x*=x2) * invfact19
    ;
  if (Math.abs(x0)>PI) { 
    x=-x;
  }

  return x;
}
class Dimension {
  float A, B, C, D;
  float F;
  float G;
  float H;
  float J;
  Dimension(float a, float b, float c, float d) {
    A=a;
    B=b;
    C=c;
    D=d;
    H =  3 * (C-B) + A-D;
    F = (B-C-C) * ( B * (B+B+C) - C * (A*3+C) )
      + D * ( A * ( A+D - 3*(B+C) ) + B * (6*B-3*C) );
    G =  9 * ( A*(C-D) + B*(D+C-B) - C*C );
    J = ( A + C - B - B ) / H;
  }
}
class Bezierr {
  Dimension[] dimension;
  Bezierr(Dimension ... xyz) {
    dimension = xyz;
  }
  private float[] points (int dim, float ... times ) {
    float[] out = new float[times.length];
    for (int i=0; i<times.length; i++) {
      float t = times[i];
      float t1 = 1.0-t;
      float A, B, C, D;
      A = dimension[dim].A;
      B = dimension[dim].B;
      C = dimension[dim].C;
      D = dimension[dim].D;
      out[i] = (A*t1 + 3*B*t)*t1*t1 + (3*C*t1 + D*t)*t*t;
    }
    return out;
  }
  boolean side(float px, float py) {
    
    /* smart decision */
    int sum=0;
    float [] points ;
    
    //x = 27*(F-H*H*time);

    //FG = x*x+4*G*G*G;
    
    boolean choice=false;
    
    if(choice){
    points = points(1, time2coords(0, px));
    for (float p : points)if ( py>p )sum++;
    }
    
    else{
    points = points(0, time2coords(1, py));
    for (float p : points)if ( px>p )sum++;
    }
    return (sum&1)!=0;
  }
  float[] time2coords(int dim, float time ) {
    float[] T = new float[3];
    float FG, x,aX, y, n, r, multn;

    float F, G, H, J;
    F = dimension[dim].F;
    G = dimension[dim].G;
    H = dimension[dim].H;
    J = dimension[dim].J;

    x = 27*(F-H*H*time);

    FG = x*x+4*G*G*G;

    aX = Math.abs(x);
    
    y = (float)Math.sqrt(Math.abs(FG));
    if (FG>=0) {
      n = (float)acbrt(Math.abs(x+y));
      multn = CBRT2*(G+G-CBRT2*n*n)/(6*H*n);
      if (x+y>=0) {
        T[0] = 0.5*multn + J;
        T[1] = T[0];
        T[2] = -multn + J;
      } else {
        T[0] = multn + J;
        T[1] = -0.5 * multn + J;
        T[2] = T[1];
      }rob=false;
    } else {
      rob=true;
      n = acbrt(aX+y);
      multn = CBRT2*(G+G-CBRT2*n*n)/(6*H*n);
      r = (fastan2(y, x)+TWO_PI)/3;
      T[0] = fastsin( -r - HALF_PI  ) * multn + J;
      T[1] = fastsin(  r + SIXTH_PI ) * multn + J;
      T[2] = fastsin( -r + SIXTH_PI ) * multn + J;
    }
    return T;
  }
}
boolean rob=false;
float x0, x1, x2, x3;
float y0, y1, y2, y3;
Bezierr bxy;
void setup() {
  size(400, 400);
  x0=300;
  x1=-50;
  x2=550;
  x3=100;
  y0=0;
  y1=133;
  y2=266;
  y3=400;

  bxy = new Bezierr(
    new Dimension(x0, x1, x2, x3), 
    new Dimension(y0, y2, y2, y3)
    );
}
void draw() {
  if (mouseX !=pmouseX || mouseY != pmouseY) {
    bxy = new Bezierr(
      new Dimension(x0, mouseX, x2, x3), 
      new Dimension(y0, mouseY, y2, y3)
      );
  }
    loadPixels();
    final int len  = pixels.length;
    for (int i=0; i<len; i++) {
      if (bxy.side(i%width, i/height)) {
        pixels[i]=#000000;
      } else { 
        pixels[i]=#ffffff;
      }
    }
    updatePixels();
    println(frameRate + " " + rob);
}
