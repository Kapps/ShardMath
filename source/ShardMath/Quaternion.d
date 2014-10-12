module ShardMath.Quaternion;
private import ShardMath.Matrix;
private import std.math;
import ShardMath.Vector;
import std.format;
import tested;

@nogc:

/// Provides a quaternion primarily used for 3D rotations. The imaginary components are excluded.
struct Quaternion  {

	/// Returns this vector as a multi-line string in the same way that an array would be.
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const {
		sink("(");
		sink.formatValue(x, fmt);
		sink(", ");
		sink.formatValue(y, fmt);
		sink(", ");
		sink.formatValue(z, fmt);
		sink(", ");
		sink.formatValue(w, fmt);
		sink(")");
	}

	@name("ToString Tests")
	unittest {
		auto quat = Quaternion(1, 2, 3, 4);
		import std.conv;
		assert(quat.text == "(1, 2, 3, 4)");
	}

@safe pure nothrow:

	/// Initializes a new Quaternion with the given components.
	this(float x, float y, float z, float w) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;
	}
	
	/// Initializes a new Quaternion with a Vector and Scalar.
	/// Params:
	/// 	vector = The elements of this Quaternion in vector form.
	/// 	scalar = The scalar value for the Quaternion.
	this(Vector3f vector, float scalar) {
		this(vector.x, vector.y, vector.z, scalar);
	}
	
	@name("Constructor Tests")
	unittest {
		Quaternion quat = Quaternion.identity;
		assert(quat.x == 0 && quat.y == 0 && quat.z == 0);
		assert(quat.w == 1);
		quat = Quaternion(1, 2, 3, 4);		
		assert(quat.x == 1 && quat.y == 2 && quat.z == 3 && quat.w == 4);
	}

	/// Calculates the magnitude of this Quaternion.
	@property float magnitude() const {
		return sqrt(magnitudeSquared);
	}	

	/// Calculates the magnitude of this Quaternion without square-rooting it.
	@property float magnitudeSquared() const {
		return (x * x) + (y * y) + (z * z) + (w * w);
	}

	/// Calculates whether this Quaternion is normalized.
	@property bool isNormalized() const {
		return abs(magnitudeSquared - 1) < 0.001f; // No need to sqrt, as we're comparing to see if it's 1.
	}

	@name("Magnitude Tests")
	unittest {
		Quaternion quat = Quaternion(1, 2, 3, 4);
		assert(quat.magnitude - sqrt(30f) < 0.01f);
		assert(quat.magnitudeSquared - 30.0f < 0.01f);
		assert(Quaternion.identity.isNormalized);
		assert(!Quaternion(1, 2, 3, 4).isNormalized);
	}

	/// Creates an identity Quaternion. When used for rotation, this indicates no rotation.
	@property static Quaternion identity() {
		return Quaternion(0, 0, 0, 1);
	}

	/// Normalizes this Quaternion, altering this instance.
	void normalizeInline() {		
		float multiplier = 1f / magnitudeSquared;
		x *= multiplier;
		y *= multiplier;
		z *= multiplier;
		w *= multiplier;
	}

	@name("Normalization Tests")
	unittest {
		Quaternion quat = Quaternion(1, 2, 3, 4);
		quat.normalizeInline();
		assert(quat.w - 0.73f < 0.01f && quat.x - 0.18f < 0.01f && quat.y - 0.36f < 0.01f && quat.z - 0.54f < 0.01f);
	}

	/// Creates a new Quaternion from the given yaw, pitch, and roll.
	/// Params:
	/// 	yaw = The angle, in radians, around the x-axis.
	/// 	pitch = The angle, in radians, around the y-axis.
	/// 	roll = The angle, in radians, around the z-axis.
	public static Quaternion fromYawPitchRoll(float yaw, float pitch, float roll) {		
		roll *= 0.5f;
		pitch *= 0.5f;
		yaw *= 0.5f;
		float rollY = sin(roll), rollX = cos(roll);
		float pitchY = sin(pitch), pitchX = cos(pitch);
		float yawY = sin(yaw), yawX = cos(yaw);
		Quaternion result;
		result.x = (yawX * pitchY * rollX) + (yawY * pitchX * rollY);
		result.y = (yawY * pitchX * rollX) - (yawX * pitchY * rollY);
		result.z = (yawX * pitchX * rollY) - (yawY * pitchY * rollX);
		result.w = (yawX * pitchX * rollX) + (yawY * pitchY * rollY);
		return result;
	}

	/// Creates a Quaternion rotating angle radians around the given axis.
	/// The axis must be a unit vector.
	/// Params:
	/// 	axis = The axis to rotate around.
	/// 	angle = The angle to rotate by.
	public static Quaternion fromAxisAngle(Vector3f axis, float angle) {
		assert(axis.isNormalized);
		float halfAngle = angle * 0.5f;
		float y = sin(halfAngle);
		float x = cos(halfAngle);
		return Quaternion(axis.x * y, axis.y * y, axis.z * y, x);
	}

	@name("Euler Tests")
	unittest {
		Quaternion quat = Quaternion.fromYawPitchRoll(1, 2, 3);
		assert(quat.w - 0.435 < 0.01f);
		assert(quat.x - 0.310 < 0.01f);
		assert(quat.y + 0.718 < 0.01f);
		assert(quat.z - 0.444 < 0.01f);
		Quaternion xQuat = Quaternion.fromAxisAngle(Vector3f(1, 0, 0), PI_2);
	}

	public Quaternion opBinary(string Op)(in Quaternion other) const if(Op == "+" || Op == "-" || Op == "/" || Op == "*") {
		static if(Op == "+" || Op == "-") {
			mixin("return Quaternion(this.x " ~ Op ~ " other.x, this.y " ~ Op ~ " other.y, this.z " ~ Op ~ " other.z, this.w " ~ Op ~ " other.w);");
		} else static if(Op == "*") {
			return Quaternion(
				(this.x * other.w) + (this.w * other.x) + (this.y * other.z) - (this.z * other.y),
				(this.y * other.w) + (this.w * other.y) + (this.z * other.x) - (this.x * other.z),
				(this.z * other.w) + (this.w * other.z) + (this.x * other.y) - (this.y * other.x),
				(this.w * other.w) - (this.x * other.x) - (this.y * other.y) - (this.z * other.z)
			);
		} else static if(Op == "/") {
			// TODO: Verify this.
			float magSquared = other.magnitudeSquared;
			float inverseMag = 1f / magSquared;
			float NX = -other.x * inverseMag, NY = -other.y * inverseMag, NZ = -other.z * inverseMag;
			float NW = other.w * inverseMag;
			return Quaternion(
				(this.x * NW) + (NX * this.w) + (this.y * NZ) - (this.z * NY),
				(this.y * NW) + (NY * this.w) + (this.z * NX) - (this.x * NZ),
				(this.z * NW) + (NZ * this.w) + (this.x * NY) - (this.y * NX),
				(this.w * NW) - ((this.x * NX) + (this.y * NY) + (this.z * NZ))
			);
		} else static assert(0, "Unsupported quaternions operator \'" ~ Op ~ "\'.");
	}

	bool opEquals(in Quaternion other) @safe const pure nothrow {
		return approxEqual(x, other.x) && approxEqual(y, other.y) && approxEqual(z, other.z) && approxEqual(w, other.w);
	}

	@name("Binary Operation Tests")
	unittest {
		Quaternion quat = Quaternion(1, 2, 3, 4);
		auto quatOrig = quat;
		quat = quat * quat;
		assert(quat.x == 8 && quat.y == 16 && quat.z == 24 && quat.w == 2);
		// TODO: Division test.
		/+Quaternion divided = quat / Quaternion(1, 2, 3, 4);
		assert(divided == Quaternion(1, 2, 3, 4));+/
		assert(quatOrig + quatOrig == Quaternion(2, 4, 6, 8));
		assert(quat - quat == Quaternion(0, 0, 0, 0));
	}

	/// Returns a normalized value of this Quaternion.
	@property Quaternion normalized() const {
		auto mag = magnitude;
		if(abs(mag) < 0.00001f || abs(mag - 1) < 0.00001f)
			return this;
		float multiplier = 1f / mag;
		return Quaternion(x * multiplier, y * multiplier, z * multiplier, w * multiplier);
	}

	@name("Normalize Tests")
	unittest {
		auto quat = Quaternion(1, 2, 3, 4);
		assert(abs(quat.normalized.magnitudeSquared - 1) < 0.001f);
	}

	/// Calculates the axis-angle components of this quaternion.
	void getAxisAngle(ref Vector3f axis, ref float angle) {
		float scale = (x * x) + (y * y) + (z * z);
		float scaleRecip = 1 /scale;
		axis = Vector3f(x * scaleRecip, y * scaleRecip, z * scaleRecip);
		angle = acos(w) * 2.0f;
	}

	/// Returns a rotation Matrix representing this Quaternion.
	Matrix4f toMatrix() const {
		Matrix4f result = void;
		float XX = this.x * this.x;
		float YY = this.y * this.y;
		float ZZ = this.z * this.z;
		float XY = this.x * this.y;
		float ZW = this.z * this.w;
		float ZX = this.z * this.x;
		float YW = this.y * this.w;
		float YZ = this.y * this.z;
		float XW = this.x * this.w;
		result.m11 = 1f - (2f * (YY + ZZ));
		result.m12 = 2f * (XY + ZW);
		result.m13 = 2f * (ZX - YW);
		result.m14 = 0f;
		result.m21 = 2f * (XY - ZW);
		result.m22 = 1f - (2f * (ZZ + XX));
		result.m23 = 2f * (YZ + XW);
		result.m24 = 0f;
		result.m31 = 2f * (ZX + YW);
		result.m32 = 2f * (YZ - XW);
		result.m33 = 1f - (2f * (YY + XX));
		result.m34 = 0f;
		result.m41 = 0f;
		result.m42 = 0f;
		result.m43 = 0f;
		result.m44 = 1f;
		return result;
	}

	@name("Conversion Tests")
	unittest {
		Quaternion quat = Quaternion(1, 2, 3, 4);
		Matrix4f mat = quat.toMatrix();
		Matrix4f expected = Matrix4f(-25, 28, -10, 0, -20, -19, 20, 0, 22, 4, -9, 0, 0, 0, 0, 1);
		assert(mat == expected);
	}

	union {
		/// The elements of this Quaternion in array form.
		float[4] elements;
		struct {
			/// The X, Y, Z, and W elements in this Quaternion.
			float x, y, z, w;
		}
	}
	
private:	
	
}