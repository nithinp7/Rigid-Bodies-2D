
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
public final class DistanceConstraint extends Constraint {

    private final float distance;
    private final Vec2 locAncA, locAncB, locAncAn, locAncBn;
    
    public DistanceConstraint(ConvPoly a, ConvPoly b, Vec2 locAncA, Vec2 locAncB, float distance) {
        super(a, b);
        this.distance = distance;
        this.locAncA = locAncA;
        this.locAncB = locAncB;
        locAncAn = locAncA.normalizeSafe();
        locAncBn = locAncB.normalizeSafe();
    }
    
    @Override
    public void render(PGraphics g) {
        Line.drawLine(a.toGlobal(locAncA), b.toGlobal(locAncB), Color.GREEN, g);
    }
    
    @Override
    public float getBias(float dt) {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA));
        
        float d = dif.mag()-distance,
              vRel = b.getVelocityLocal(locAncB).sub(a.getVelocityLocal(locAncA))
                        .dot(dif.normalize());
        
        return d/dt;// - vRel;//-0.01f*d/dt - vRel;
    }
    
    @Override
    public float[] getJacobian() {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA)),
             difN = dif.normalizeSafe();
        return new float[] {
            -difN.x, -difN.y, -a.orient(locAncAn).cross(difN),
            difN.x, difN.y, b.orient(locAncBn).cross(difN)
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
