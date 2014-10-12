module ShardMath.BoundingFrustum;
import ShardMath.Matrix;
import ShardMath.Vector;
import ShardMath.BoundingBox;
import std.typecons;

/*/// Represents a frustum in 3D space, generally used for view culling.
@disable class BoundingFrustum  {

public:	

	/// Initializes a new instance of the BoundingFrustum object.
	/// Params:
	///		Matrix = The world matrix multiplied by the view matrix multiplied by the projection matrix.
	this(Matrix4f WorldViewProj) {
		this._Matrix = WorldViewProj;
	}

	/// Determines whether this BoundingFrustum contains the given object.
	/// Params:
	/// 	Object = The object to check for containment.
	ContainmentType Overlap(BoundingBox Object) {
		throw new NotImplementedError("Overlap");
	}

	/// Ditto
	ContainmentType Overlap(Vector3f Point) {
		throw new NotImplementedError("Overlap");
	}

	/// Gets or sets the matrix used for this BoundingFrustum. This is generally the view matrix multiplied by the projection matrix.
	@property Matrix4f Matrix() const {
		return _Matrix;
	}

	/// Ditto
	@property void Matrix(Matrix4f Value) {
		this._Matrix = Value;
	}
	
private:
	Matrix4f _Matrix;
}*/