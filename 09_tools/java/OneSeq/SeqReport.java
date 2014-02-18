/**
 *
 * @author Ivica 2005
 */

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.image.*;
import javax.imageio.*;
import javax.swing.ImageIcon;


public class SeqReport {

    final float MAXCOLOR = 1f;  // for Ivana color scheme
    float[] x = new float[10];
    int[] y = new int[10];
    String[] residue = new String[10];
    float minX = 1000.f;
    float maxX = -1000.f;
    private static final int TT_NUMBER = StreamTokenizer.TT_NUMBER;
    private static final int TT_EOF = StreamTokenizer.TT_EOF;
    private static final int TT_WORD = StreamTokenizer.TT_WORD;
    private boolean missing = false;
    public void readPoints(File file, int start, int end) {

       	int i = 0;
	int cnt = 0;
	int arrayIndex = 0;
	float tempInput = 0;
	int tempExtra = 0;
	String line;
	String stringInput = "";
	StreamTokenizer tok;

	try {

	    Reader reader = new BufferedReader( new FileReader(file));
	    tok = new StreamTokenizer(reader);
	    tok.eolIsSignificant(true);
	    skipLeadingNewLines(tok);

	    while (tok.ttype != TT_EOF) {
		if (tok.nextToken() != TT_NUMBER)
		    throw new Exception ("Coverage missing at input line " + tok.lineno());
		tempInput = (float)tok.nval;
	      	if (tok.nextToken() != TT_WORD)
		    throw new Exception ("Residue missing at input line " + tok.lineno());
		stringInput = tok.sval;
		if (tok.nextToken() != TT_NUMBER)
		    throw new Exception ("Coverage missing at input line " + tok.lineno());
		tempExtra = (int)tok.nval;

		if (cnt >= start && cnt <= end) {

		    if (x.length == arrayIndex) {
			//dimension of ranks array must be equal to input file size
			float[] temp = new float[arrayIndex+1];
			System.arraycopy(x, 0, temp,0, arrayIndex);
			x = temp;
			int[] tempY = new int[arrayIndex+1];
			System.arraycopy(y, 0, tempY,0, arrayIndex);
			y = tempY;
			String[] strTemp = new String[arrayIndex+1];
			System.arraycopy(residue, 0, strTemp,0, arrayIndex);
			residue = strTemp;

		    }

		    x[arrayIndex] = tempInput;
		    residue[arrayIndex] = stringInput;
		    y[arrayIndex] = tempExtra;
		    arrayIndex++;
		}
		cnt++;

		skipLeadingNewLines(tok);
	    }

	    int countMissing = 0;
	    for (int j = 0; j < y.length - 1; j++) {
	    	if ( (int)(y[j + 1] - y[j]) != 1) {
		    missing = true;
		    countMissing ++;
	    	}
	    }
	    if (missing) {
	    	float[] modX = new float[x.length + 3*countMissing];
	    	int[] modY = new int[x.length + 3*countMissing];
	    	String[] modResidue = new String[x.length + 3*countMissing];
	    	int k = 0;
	    	int j = 0;

		for (j = 0; j < y.length - 1; j++) {
		    modX[k] = x[j];
		    modResidue[k] =	residue[j];
		    modY[k] = y[j];
		    k++;
		    for (int iii = 0; iii < 3; iii++) {

			if ( (int)(y[j + 1] - y[j]) != 1) {
			    modX[k] = -10;
			    modResidue[k] = "dots";
			    modY[k] = -10;
			    k++;
			}
		    }

		}
		modX[k] = x[j];
		modResidue[k] =	residue[j];
		modY[k] = y[j];
		x = modX;
		residue = modResidue;
		y = modY;
	    }


	}
	catch (Exception e) {
	    System.out.println("Exception " + e.getMessage() + " in readRanks()");
	}


    }

    public static void skipLeadingNewLines(StreamTokenizer t) throws IOException {
	while (t.nextToken() == t.TT_EOL);
	t.pushBack();
    }

    class MapPanel extends JComponent {
	//this class does all the drawing

	float width = 500;//0;  // frame size
	float height = 500;//0;
	float margin =25;//20;
	float marginTop = 25;//40;
	float squareX = 20;  // legend bar width

	String title;
	Font numberFont = new Font("Monospaced",Font.PLAIN,20);
	Font residueFont = new Font("Monospaced",Font.PLAIN,10);


	MapPanel(String title) {
	    this.title = title;
	    setPreferredSize(new Dimension((int)width, (int)height));
	}

	MapPanel(String title, float width, float height) {
	    this.title = title;
	    this.width = width;
	    this.height = height;
	    setPreferredSize(new Dimension((int)width, (int)height));
	}

	public Color setIvanaColor(double cvg) {
	    //color scheme used by Ivana
	    double red, blue, green, ratio, x;
	    int  bin_size;
	    int range = 20;
	    int N = 5;
	    int color_index;
	    red = green = blue = 254;

	    if (cvg < 0.05) {
		return new Color(254,(int)(0.83*254), (int)(0.17*254));
	    } 	else {

		color_index = (int) (range * cvg);
		bin_size = (int)Math.round ( (double)(range - 1) / N);

		if (cvg <= 0.25) {
		    ratio =  (double)(bin_size - color_index + 1)/bin_size;
		    if(ratio > 1)
			ratio = 1;
		    red = ratio * 254;
		    green = blue = 0.;

		} else {
		    ratio = ( color_index  - range / N) / ( (double) (N-1.)/ N * range);
		    if(ratio > 1)
			ratio = 1;
		    red = ratio * 254;
		    green = blue = red;
		}
	    }
	    int re, bl, gr;
	    re = (int)red;
	    bl = (int)blue;
	    gr = (int)green;
	    return new Color((int)red, (int)blue, (int)green);

	}

	//coloring of data points
	public Color myColor(double norm) {
	    if (norm < 0.95)     //to avoid starting end ending in red
		return Color.getHSBColor((float)norm, 1, 1);
	    else
		return Color.getHSBColor(0.95f, 1, 1);
	}

	// to draw the color scheme
	public void drawLegend(Graphics2D g, float x, float y) {

	    int i;
	    float hue;
	    Color color;
	    Font defaultFont = g.getFont();
	    FontMetrics fm = g.getFontMetrics(numberFont);
	    String top, bottom;

	    g.setFont(numberFont);

	    for (i = 0; i < (int)y; i++) {

		//draw color bar as a series of tiny rectangles
		hue = MAXCOLOR *i / y;
		if (hue < 0.)
		    color = Color.black;
		else
		    color = setIvanaColor(hue);
		//	color = myColor(hue);
		g.setPaint(color);
		Rectangle2D.Float rectangle = new Rectangle2D.Float(x, i,squareX,1);
		g.fill(rectangle);
	    }
	    g.setPaint(Color.black);
	    g.drawString("1",x+squareX+5,i + 6);
	    i=0;

	    g.drawString(Integer.toString(i),x+squareX+5,i + 6);

	    g.setPaint(Color.black);
	    g.setFont(defaultFont);

	}


	public void paintComponent(Graphics _g) {
	    drawSequence(_g,getWidth(), getHeight(),0);
	}

	public void drawSequence(Graphics _g, float width, float height, int shift) {
	    //put all this together and paint
	    int i, j, titleCenter;
	    int position = 0;
	    int max = x.length/50 + 1;
	    float delta, stripe, lineWidth;
	    float cov=1.1f;//0;
	    float spacing;

	    Font titleFont = new Font("Monospaced",Font.PLAIN,30);

	    FontMetrics fm = _g.getFontMetrics(titleFont);
	    Rectangle2D.Float rectangle;
	    GeneralPath boundary = new GeneralPath();
	    Graphics2D g = (Graphics2D)_g;

	    g.setColor(Color.black);
	    g.setFont(titleFont);
	    titleCenter = fm.stringWidth(title);
	    stripe = (width - 2 * margin) / 50;
	    spacing = (height - margin - marginTop) / 19;
	    spacing = 23.158f;
	    //getting nice lines
	    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
			       RenderingHints.VALUE_ANTIALIAS_ON);
	    g.setStroke(new BasicStroke(1.0f));

	    delta = 0;int once = 0;
	    int countResidues = 0;
	    for (i = 0; i < x.length; i++) {
		if (i%50 == 0)
		    delta= 0;
		g.setPaint(setIvanaColor(x[i]));
		if(x[i] / MAXCOLOR > cov)
		    g.setPaint(Color.blue);
		if (i != x.length-1 ||  ((i+1)%50 == 0) )
		    lineWidth = stripe*1.1f;//1.25f; //stripe width slightly enlarged
		//to avoid white stripes between adjacent stripes
		else
		    lineWidth = stripe;
		if ( (i != x.length-1) && x[i+1] < 0)
		    lineWidth = stripe;

		rectangle = new Rectangle2D.Float(margin + delta, marginTop + position,lineWidth,spacing);


		if(x[i] >= 0) {
		    g.fill(rectangle);
		    g.setPaint(Color.black);
		    boundary.moveTo(margin + delta, marginTop + position);
		    boundary.lineTo(margin + delta + lineWidth, marginTop + position);
		    boundary.moveTo(margin + delta, marginTop + position + spacing);
		    boundary.lineTo(margin + delta + lineWidth, marginTop + position + spacing);

		    g.setPaint(Color.black);

		    g.setFont(residueFont);
		    String axes = "";
		    axes = residue[i];
		    //	if (i < x.length - 5)
		    g.drawString(axes,margin + delta + 1.5f, marginTop + position + stripe + 25) ;
		    if(countResidues%10 == 0){ 
			axes = Integer.toString(y[i]);
			g.drawString(axes,margin + delta, marginTop + position + stripe + 33.5f) ;
		    }
		    countResidues++;
		}
		else {
		    Ellipse2D circle = new Ellipse2D.Float(margin + delta + lineWidth / 2 - 1.5f,
							   marginTop + position + spacing/2, 3, 3);
		    g.setPaint(Color.black);
		    g.fill(circle);


		}

		delta +=stripe;



		if ((i+1)%50 == 0) {
		    position += spacing + 32; // seq width + space between sequences
		}

	    }


	    g.draw(boundary);
	    g.translate(width - 2*margin, marginTop);

	}
    }   //end of class MapPanel


    public static void main(String[] args) throws Exception {
	if (args.length < 2) {
	    System.out.println("Usage: java SeqReport <coverage file> \"figure name\" <initial number>  <end number>");
	    System.out.print("<coverage file> is a 3 column file ");
	    System.out.println("of form \"coverage residue_name residue_number\", ordered according to residue position on sequence");
	    System.out.println("If no numbers given, the whole coverage file will be used");
	    System.out.println("If <end number> not given, EOF is assumed to be the end");
	    System.exit(0);
	}

	int w, h, shift, end;
	float dimX, dimY;
    	Graphics2D g;
        BufferedImage image;

	// create a new Plot object
	SeqReport plot = new SeqReport();

	//read input file (should have 3 columns: x,y and z
	if (args.length >= 3 && Integer.parseInt(args[2]) > 0) {
	    shift = Integer.parseInt(args[2]) - 1;
	}
	else
	    shift = 0;
	if (args.length == 4 && Integer.parseInt(args[3]) > Integer.parseInt(args[2]))
	    end = Integer.parseInt(args[3]) - 1;
	else
	    end = 10000;

	plot.readPoints(new File(args[0]), shift, end);
	//panel used for drawing
	//   MapPanel map = plot.new MapPanel(args[1]);

        dimY = 2*(plot.x.length/50 + 1) - 1;
	dimY = 3*(plot.x.length/50 + 1) - 1;
	dimY *= 23.157f;
     
	dimY = 1.7f * (plot.x.length/50 + 1) * 32  + 10;
	MapPanel map = plot.new MapPanel(args[1], 500, dimY);

        w = (int)map.width;
        h = (int)map.height + 25;
        // 25 added for png images production
        // remove 25 im eps is going to be needed

       	image = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        g = image.createGraphics();
       	g.setColor(Color.white);
       	g.fillRect(0, 0, w, h);
      	map.drawSequence(g,w,h,shift);
       	ImageIO.write(image, "png", new File(args[1]+".png"));


    }

}   //end of class Plot

class TreeFrame extends JFrame {
    // JFrame which draws JComponent
    TreeFrame(JComponent cmp) {
    	getContentPane().setBackground(Color.white);
	getContentPane().add(cmp);
	pack();
	addWindowListener(new WindowAdapter(){
    		public void windowClosing(WindowEvent evt){
		    System.exit(0);
      		}
	    });
    	setVisible(true);
    }
}
