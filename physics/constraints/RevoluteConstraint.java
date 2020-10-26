
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
public final class RevoluteConstraint extends Constraint {

    private final Vec2 locAncA, locAncB;
    
    public RevoluteConstraint(ConvPoly a, ConvPoly b, Vec2 locAncA, Vec2 locAncB) {
        super(a, b);
        this.locAncA = locAncA;
        this.locAncB = locAncB;
    }
    
    @Override
    public void render(PGraphics g) {
        Line.drawLine(a.toGlobal(locAncA), b.toGlobal(locAncB), Color.GREEN, g);
    }
    
    @Override
    public float getBias(float dt) {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA));
        
        float d = dif.mag(),
              vRel = b.getVelocityLocal(locAncB).sub(a.getVelocityLocal(locAncA))
                        .dot(dif.normalize());
        
        return d/dt;//getRestitution()*vRel - 0.2f*d/dt;
    }
    
    @Override
    public float[] getJacobian() {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA));
        if(dif.equals(Vec2.ZERO)) 
            return null;
        Vec2 difN = dif.normalize();
        return new float[] {
            -difN.x, -difN.y, -a.orient(locAncA).cross(difN),
            difN.x, difN.y, b.orient(locAncB).cross(difN)
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
