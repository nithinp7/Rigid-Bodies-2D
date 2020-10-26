
package physics.constraints;

import geometry.visual.Line;
import java.awt.Color;
import physics.shapes.ConvPoly;
import static processing.core.PApplet.*;
import processing.core.PGraphics;

/**
 *
 * @author Nithin
 */
public final class GearConstraint extends Constraint {

    //TODO: implement gear ratio
    
    //angleOffset is the target angle of b relative to a
    private final float angleOffset;
    
    public GearConstraint(ConvPoly a, ConvPoly b, float angleOffset) {
        super(a, b);
        this.angleOffset = angleOffset%(2*PI);
    }
    
    public GearConstraint(ConvPoly a, ConvPoly b) {
        this(a, b, b.getAngle()-a.getAngle());
    }
    
    @Override
    public void render(PGraphics g) {
        Line.drawLine(a.getLocation(), b.getLocation(), Color.BLUE, g);
    }
    
    @Override
    public float getBias(float dt) {
        return -0.5f*(b.getAngle() - a.getAngle() - angleOffset)/dt;
    }
    
    @Override
    public float[] getJacobian() {
        return new float[] {
           0, 0, 1,
           0, 0, -1
        };
    }
    
    @Override
    public float getLo() {
        return Float.NEGATIVE_INFINITY;
    }
    
    @Override
    public float getHi() {
        return Float.POSITIVE_INFINITY;
    }
}
