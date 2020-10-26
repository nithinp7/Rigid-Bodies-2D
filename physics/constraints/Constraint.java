
package physics.constraints;

import physics.shapes.ConvPoly;
import static processing.core.PApplet.min;
import processing.core.PGraphics;

/**
 *
 * @author Nithin
 */
public abstract class Constraint {
    
    public final ConvPoly a, b;
    
    public Constraint(ConvPoly a, ConvPoly b) {
        this.a=a;
        this.b=b;
    }
    
    public abstract void render(PGraphics g);
    
    public float getRestitution() {
        return min(a.e, b.e);
    }
    
    public abstract float getBias(float dt);
    
    public abstract float[] getJacobian();
    
    public abstract float getLo();
    
    public abstract float getHi();
}
