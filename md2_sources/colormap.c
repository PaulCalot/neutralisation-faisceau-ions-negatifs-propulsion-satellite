/* colormap.c	(c) 1997 Cam Abrams
 * University of California,  Berkeley
 * Dept. of Chemical Engineering
 * Graves Group
 */

#include "colormap.h"

float hue_max_ = 1.0;
float hue_min_ = 0.0;


/*#include "cam_colors.cf"*/
color_t black = {0.0, 0.0, 0.0};
color_t white = {1.0, 1.0, 1.0};
color_t red = {1.0, 0.0, 0.0};
color_t orange = {1.0, 0.5, 0.0};
color_t green = {0.0, 1.0, 0.0};
color_t blue = {0.0, 0.0, 1.0};
color_t darkgold1 = { 0.933, 0.706, 0.059};
color_t peagreen = { 0.2, 1.0, 0.2 };
color_t darkpeagreen = { 0.1, 0.5, 0.1 };
color_t yellow = { 1.0, 1.0, 0.0 };
color_t violet = { 1.0, 0.0, 1.0 };
color_t violet2 = { 0.85, 0.3, 0.85 };
color_t bronze = { 0.545, 0.396, 0.031 };
color_t grey1 = { 0.5, 0.5, 0.5 };
color_t grey2 = { 0.25, 0.25, 0.25 };
color_t grey3 = { 0.125, 0.125, 0.125 };
color_t cyan = { 0.0, 1.0, 1.0 };
color_t SlateBlue = { 0.4157, 0.3529, 0.8039 };
color_t navy = { 0, 0, 0.50196 };
color_t khaki = { 0.94118, 0.90196, 0.54902 };
color_t grey4 = { 0.07843, 0.07843, 0.07843 };

hue_t color_maxHue ( color_t *c )
{
    if (c->r >= c->b && c->r >= c->g) return C_R;
    if (c->g >= c->b && c->g >= c->r) return C_G;
    if (c->b >= c->r && c->b >= c->g) return C_B;
    return C_R;
}

hue_t color_minHue ( color_t *c )
{
    if (c->r <= c->b && c->r <= c->g) return C_R;
    if (c->g <= c->b && c->g <= c->r) return C_G;
    if (c->b <= c->r && c->b <= c->g) return C_B;
    return C_R;
}

float color_maxVal ( color_t *c )
{
    return color_hueVal(c, color_maxHue(c));
}

float color_minVal ( color_t *c )
{
    return color_hueVal(c, color_minHue(c));
}

color_t *color_init ( color_t *res, float r, float g, float b )
{
    if (!res) return NULL;
    res->r = r;
    res->g = g;
    res->b = b;
    return res;
}

color_t *colorPtr (char * colorStr)
{
    if (!colorStr) return NULL;
/*    #include "cam_str2clr.cf"*/
    if (!strcmp(colorStr, "black") || !strcmp(colorStr, "BLACK")) return &black;
    if (!strcmp(colorStr, "white") || !strcmp(colorStr, "WHITE")) return &white;
    if (!strcmp(colorStr, "red") || !strcmp(colorStr, "RED")) return &red;
    if (!strcmp(colorStr, "orange") || !strcmp(colorStr, "ORANGE")) return &orange;
    if (!strcmp(colorStr, "green") || !strcmp(colorStr, "GREEN")) return &green;
    if (!strcmp(colorStr, "blue") || !strcmp(colorStr, "BLUE")) return &blue;
    if (!strcmp(colorStr, "darkgold1") || !strcmp(colorStr, "DARKGOLD1")) return &darkgold1;
    if (!strcmp(colorStr, "peagreen") || !strcmp(colorStr, "PEAGREEN")) return &peagreen;
    if (!strcmp(colorStr, "darkpeagreen") || !strcmp(colorStr, "DARKPEAGREEN")) return &darkpeagreen;
    if (!strcmp(colorStr, "yellow") || !strcmp(colorStr, "YELLOW")) return &yellow;
    if (!strcmp(colorStr, "violet") || !strcmp(colorStr, "VIOLET")) return &violet;
    if (!strcmp(colorStr, "violet2") || !strcmp(colorStr, "VIOLET2")) return &violet2;
    if (!strcmp(colorStr, "bronze") || !strcmp(colorStr, "BRONZE")) return &bronze;
    if (!strcmp(colorStr, "grey1") || !strcmp(colorStr, "GREY1")) return &grey1;
    if (!strcmp(colorStr, "grey2") || !strcmp(colorStr, "GREY2")) return &grey2;
    if (!strcmp(colorStr, "grey3") || !strcmp(colorStr, "GREY3")) return &grey3;
    if (!strcmp(colorStr, "cyan") || !strcmp(colorStr,"CYAN")) return &cyan;
    if (!strcmp(colorStr, "SlateBlue") || !strcmp(colorStr,"SLATEBLUE")) return &SlateBlue;
    if (!strcmp(colorStr, "navy") || !strcmp(colorStr,"NAVY")) return &navy;
    if (!strcmp(colorStr, "khaki") || !strcmp(colorStr,"KHAKI")) return &khaki;
    if (!strcmp(colorStr, "grey4") || !strcmp(colorStr,"GREY4")) return &grey4;

    return NULL;
}

char tmp_colorStr[255];
char * colorStr (color_t * c)
{
    char * rv = tmp_colorStr;
    rv[0] = '\0';
/*    #include "cam_clr2str.cf"*/
    if (c == &black) strcpy(rv, "black");
    if (c == &white) strcpy(rv, "white");
    if (c == &red) strcpy(rv, "red");
    if (c == &orange) strcpy(rv, "orange");
    if (c == &green) strcpy(rv, "green");
    if (c == &blue) strcpy(rv, "blue");
    if (c == &darkgold1) strcpy(rv, "darkgold1");
    if (c == &peagreen) strcpy(rv, "peagreen");
    if (c == &darkpeagreen) strcpy(rv, "darkpeagreen");
    if (c == &yellow) strcpy(rv, "yellow");
    if (c == &violet) strcpy(rv, "violet");
    if (c == &violet2) strcpy(rv, "violet2");
    if (c == &bronze) strcpy(rv, "bronze");
    if (c == &grey1) strcpy(rv, "grey1");
    if (c == &grey2) strcpy(rv, "grey2");
    if (c == &grey3) strcpy(rv, "grey3");
    if (c == &cyan) strcpy(rv, "cyan");
    if (c == &SlateBlue) strcpy(rv, "SlateBlue");
    if (c == &navy) strcpy(rv, "navy");
    if (c == &khaki) strcpy(rv, "khaki");
    if (c == &grey4) strcpy(rv, "grey4");

    return rv;
}

hue_t color_nextHue ( hue_t Hue, int flow)
{
    if (flow != 1 && flow != -1) return Hue;
    if (Hue != C_R && Hue != C_B && Hue != C_G) return Hue;
    if (Hue == C_R && flow == 1) return C_G;
    if (Hue == C_R && flow == -1) return C_B;
    if (Hue == C_G && flow == 1) return C_B;
    if (Hue == C_G && flow == -1) return C_R;
    if (Hue == C_B && flow == 1) return C_R;
    if (Hue == C_B && flow == -1) return C_G;
    return C_G;
}

float color_hueVal ( color_t *c,  hue_t Hue )
{
    if (Hue != C_R && Hue != C_B && Hue != C_G) return 0.0;
    if (Hue == C_R) return c->r;
    if (Hue == C_G) return c->g;
    if (Hue == C_B) return c->b;
    return C_G;
}

float color_color2x ( color_t *c )
{
    if (color_maxHue(c) == C_R)
    {
	float f = c->g - c->b;
	if (f < 0.0) return (1.0 + ONE_SIXTH * f) * (hue_max_ - hue_min_);
	return ONE_SIXTH * f * (hue_max_ - hue_min_);  
    }
    if (color_maxHue(c) == C_G)
	return (ONE_THIRD + ONE_SIXTH * (c->b - c->r)) * (hue_max_ - hue_min_);
  
    if (color_maxHue(c) == C_B) 
	return (TWO_THIRDS + ONE_SIXTH * (c->r - c->g)) * (hue_max_ - hue_min_);

    return 0.0;  
}

color_t *color_x2color (color_t *c,  float x)
{
    double i,  *ip = &i;
    color_init(c, 0.0, 0.0, 0.0);
    x = modf(x, ip);
    if (x < 0.0) x += 1.0;
    if (x > (1.0 - ONE_SIXTH)) x -= 1.0;
    if (x <= ONE_SIXTH || x >= (1.0 - ONE_SIXTH)) {
	c->r = hue_max_;
	if (x < 0.0) {
	    c->g = hue_min_;
	    c->b = -6*x*(hue_max_ - hue_min_);
	}
	else {
	    c->b = hue_min_;
	    c->g = 6*x*(hue_max_ - hue_min_);
	}
    }
    else if (x <= 0.5) {
	c->g = hue_max_;
	if (x < ONE_THIRD) {
	    c->b = hue_min_;
	    c->r = 6*(ONE_THIRD - x)*(hue_max_ - hue_min_);
	}
	else {
	    c->r = hue_min_;
	    c->b = 6*(x - ONE_THIRD)*(hue_max_ - hue_min_);
	}
    }
    else {
	c->b = hue_max_;
	if (x < TWO_THIRDS) {
	    c->r = hue_min_;
	    c->g = 6*(TWO_THIRDS - x)*(hue_max_ - hue_min_);
	}
	else {
	    c->g = hue_min_;
	    c->r = 6*(x - TWO_THIRDS)*(hue_max_ - hue_min_);
	}
    }
    return c;
}

void color_printf (color_t *c) {
    if (!c) printf("0.0 0.0 0.0");
    printf("%.3f %.3f %.3f", c->r, c->g, c->b);
}

color_t *color_scalMultColor (color_t *c, double x)
{
    return color_init(c, c->r*x, c->g*x, c->b*x);
}

color_t *color_mapFloat ( float f0, color_t *c0, float f1, color_t *c1, 
			 float f, color_t *c, int flow)
{
    float M1, m1, M0, m0;
    float m = (f1 - f)/(f1 - f0);
    float x0 = color_color2x(c0);
    float x1 = color_color2x(c1);
    float dx = x1 - x0;
    float xc = 0.0;
    
    M1 = color_maxVal(c1);
    M0 = color_maxVal(c0);
    m1 = color_minVal(c1);
    m0 = color_minVal(c0);
    
    /* check constraints on c0 and c1 */
    if ((M1 != M0) || (m1 != m0)) return color_init(c, 0.0, 0.0, 0.0);
    
    /* if the float is not within the range (f0,f1) then assign the slope
     * so that if f < f0, xc = x0, and if f > f1, xc = x1.
     */
    
    if (f > f0 && f < f1)
    {
	if (flow == 1) xc = m*dx + x0;
	else xc = x1 - m*dx;
    }
    else if (f <= f0)
	xc = x0;
    else if (f >= f1)
	xc = x1;
    
    
    return color_x2color (c, xc);

}
