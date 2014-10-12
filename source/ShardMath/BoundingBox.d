module ShardMath.BoundingBox;
import ShardMath.Vector;
import tested;

@nogc:

/// Determines the type of containment between two objects, with values being arranged by most containment (fully), to least containment (disjoint).
enum OverlapType {
	///
	disjoint,
	///
	intersects,
	///
	contains
}

/// Represents a 3D Axis-Aligned Bounding Box (AABB).
struct BoundingBox  {

@safe pure nothrow:

	/// The minimum coorindates for this box.
	Vector3f min;
	/// The maximum coordinates for this box.
	Vector3f max;

	/// Creates a BoundingBox with the given coordinates.
	/// Params:
	/// 	min = The minimum coordinates for this BoundingBox.
	/// 	max = The maximum coordinates for this BoundingBox.
	this(Vector3f min, Vector3f max) {
		this.min = min;
		this.max = max;
		assert(min.x <= max.x && min.y <= max.y && min.z <= max.z);
	}
	
	/// Checks whether the two BoundingBoxes collide with each other.
	/// Params:
	/// 	other = The other BoundingBox to check for collision.
	bool intersects(BoundingBox other) const {
		if(max.x < other.min.x || min.x > other.max.x)
			return false;
		if(max.y < other.min.y || min.y > other.max.y)
			return false;
		return max.z >= other.min.z && min.z <= other.max.z;
	}

	/// Determines how the two BoundingBoxes are contained within each other.
	/// Params:
	/// 	other = The other BoundingBox to check for collision.
	OverlapType overlap(BoundingBox other) const {
		if(max.x < other.min.x || min.x > other.max.x)
			return OverlapType.disjoint;
		if(max.y < other.min.y || min.y > other.max.y)
			return OverlapType.disjoint;
		if(max.z < other.min.z || min.z > other.max.z)
			return OverlapType.disjoint;
		if(min.x <= other.min.x && max.x >= other.max.x && min.y <= other.min.y && max.y >= other.max.y && min.z <= other.min.z && max.z >= other.max.z)
			return OverlapType.contains;
		return OverlapType.intersects;
	}

	bool opEquals(const ref BoundingBox other) const {
		return min == other.min && max == other.max;
	}

	/// Returns the difference between the maximum and minimum boundaries of this BoundingBox.
	@property Vector3f diff() const {
		return max - min;
	}

	@name("BoundingBox Tests")
	unittest {
		BoundingBox first = BoundingBox(Vector3f(0.33, 0.33, 0.33), Vector3f(1, 1, 1));
		BoundingBox firstEqual = first;
		BoundingBox second = BoundingBox(Vector3f(0, 0, 0), Vector3f(2, 2, 2));
		BoundingBox third = BoundingBox(Vector3f(0, 0, 0), Vector3f(1, 1, 1));

		//assert(first == firstEqual);
		assert(first.intersects(second));
		assert(second.intersects(first));
		assert(second.overlap(first) == OverlapType.contains);
		assert(second.overlap(third) == OverlapType.contains);
		assert(third.overlap(second) == OverlapType.intersects);
		assert(first.overlap(second) == OverlapType.intersects);
	}
	
}