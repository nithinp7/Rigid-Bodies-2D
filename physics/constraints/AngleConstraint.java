
package physics.constraints;

import geometry.visual.Line;
import java.awt.Color;
import math.Vec2;
import physics.shapes.ConvPoly;
import processing.core.PGraphics;

/**
 *
 * @author Nithin
 */
public class AngleConstraint extends Constraint {
    
    //angleOffset is the target angle of b relative to a
    private final Vec2 locAncA, locAncB, locNormA, locNormB;
    
    public AngleConstraint(ConvPoly a, ConvPoly b,
                           Vec2 locAncA, Vec2 locAncB,
                           Vec2 globalNormal) {
        
        super(a, b);
        
        this.locAncA = locAncA;
        this.locAncB = locAncB;
        
        globalNormal = globalNormal.normalize();
        
        locNormA = a.disorient(globalNormal);
        locNormB = b.disorient(globalNormal);
    }
    
    /*
    public AngleConstraint(ConvPoly a, ConvPoly b) {
        this(a, b, b.getAngle()-a.getAngle());
    }
    */
    @Override
    public void render(PGraphics g) {
        Line.drawLine(a.getLocation(), b.getLocation(), Color.BLUE, g);
    }
    
    @Override
    public float getBias(float dt) {
        return -0.5f*(b.getAngle() - a.getAngle() - 0/*angleOffset*/)/dt;
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
