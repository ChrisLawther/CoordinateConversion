import XCTest
import CoreLocation
@testable import CoordinateConversion

final class CoordinateConverstionTests: XCTestCase {
    func testCatAndFiddle() throws {

        // Cat and Fiddle pub in Derbyshire
        let (e, n) = (400100.0, 371875.0)
        let cat_n_fiddle = CLLocation(easting: e, northing: n)

        let bng = cat_n_fiddle.toBNG()

        // Tiny errors after round-trip, hopefully
        XCTAssertEqual(e, bng.eastings, accuracy: 0.01)
        XCTAssertEqual(n, bng.northings, accuracy: 0.01)

    }
}
