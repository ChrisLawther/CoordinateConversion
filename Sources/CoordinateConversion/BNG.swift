import CoreLocation

/// Describes a location using the British National Grid (aka Ordnance Survey National Grid)
/// system
///
/// https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid
public struct BNGLocation {
    public let easting: Double
    public let northing: Double

    public init(easting: Double, northing: Double) {
        self.easting = easting
        self.northing = northing
    }

    /// Initialise a British National Grid location from a CLLocation
    /// - Parameter location: The `CLLocation`
    public init(_ location: CLLocation) {
        let coordinate = location.coordinate

        let lat1 = coordinate.latitude * Double.pi / 180.0
        let lon1 = coordinate.longitude * Double.pi / 180.0

        let (a1, b1) = (6378137.000, 6356752.3141)

        var e2 = 1 - (b1 * b1) / (a1 * a1)
        let nu1 = a1 / sqrt(1 - e2 * pow(sin(lat1),2))

        // First convert to cartesian from spherical polar coordinates
        var H = 0.0               // Third spherical coord
        let x1 = (nu1 + H) * cos(lat1) * cos(lon1)
        let y1 = (nu1 + H) * cos(lat1) * sin(lon1)
        let z1 = ((1 - e2) * nu1 + H) * sin(lat1)

        // Perform Helmut transform (to go between GRS80 (_1) and Airy 1830 (_2))
        let s = 20.4894E-6      // The scale factor -1

        // The translations along x,y,z axes respectively
        let (tx, ty, tz) = (-446.448, 125.157, -542.060)

        // The rotations along x,y,z respectively, in seconds
        let (rxs, rys, rzs) = (-0.1502, -0.2470, -0.8421)

        // ... in radians
        let (rx, ry, rz) = (rxs * Double.pi / (180 * 3600.0), rys * Double.pi / (180 * 3600.0), rzs * Double.pi / (180 * 3600.0))

        let x2 = tx + (1 + s) * x1 + (-rz) * y1 + (ry) * z1
        let y2 = ty + (rz) * x1 + (1 + s) * y1 + (-rx) * z1
        let z2 = tz + (-ry) * x1 + (rx) * y1 + (1 + s) * z1

        // Back to spherical polar coordinates from cartesian
        // Need some of the characteristics of the new ellipsoid
        // The GSR80 semi-major and semi-minor axes used for WGS84(m)
        let (a, b) = (6377563.396, 6356256.909)

        e2 = 1 - (b * b) / (a * a)              // The eccentricity of the Airy 1830 ellipsoid
        let p = sqrt(x2 * x2 + y2 * y2)

        // Lat is obtained by an iterative proceedure:
        var lat = atan2(z2, (p * (1 - e2)))     // Initial value
        var latold = 2 * Double.pi

        var nu = 0.0
        while abs(lat - latold) > 1e-16 {
            (lat, latold) = (latold, lat)
            nu = a / sqrt(1 - e2 * pow(sin(latold),2))
            lat = atan2(z2 + e2 * nu * sin(latold), p)
        }

        // Lon and height are then pretty easy
        let lon = atan2(y2, x2)
        H = p / cos(lat) - nu

        // E, N are the British national grid coordinates - eastings and northings
        let F0 = 0.9996012717                   // scale factor on the central meridian
        let lat0 = 49 * Double.pi / 180.0            // Latitude of true origin (radians)
        let lon0 = -2 * Double.pi / 180.0            // Longtitude of true origin and central meridian (radians)
        let (N0, E0) = (-100000.0, 400000.0)    // Northing & easting of true origin (m)
        let n = (a - b) / (a + b)

        // meridional radius of curvature
        let rho = a * F0 * (1 - e2) * pow((1 - e2 * pow(sin(lat),2)),-1.5)
        let eta2 = nu * F0 / rho - 1

        let M1 = (1 + n + (5.0 / 4.0) * pow(n,2) + (5.0 / 4.0) * pow(n,3)) * (lat-lat0)
        let M2 = (3.0 * n + 3.0 * pow(n,2) + (21.0 / 8.0) * pow(n,3)) * sin(lat - lat0) * cos(lat + lat0)
        let M3 = ((15.0 / 8.0) * pow(n,2) + (15.0 / 8.0) * pow(n,3)) * sin(2.0 * (lat-lat0)) * cos(2.0 * (lat+lat0))
        let M4 = (35.0 / 24.0) * pow(n,3) * sin(3.0 * (lat-lat0)) * cos(3.0 * (lat+lat0))

        // meridional arc
        let M = b * F0 * (M1 - M2 + M3 - M4)

        let I = M + N0
        let II = nu * F0 * sin(lat) * cos(lat) / 2.0
        let III = nu * F0 * sin(lat) * pow(cos(lat),3) * (5 - pow(tan(lat),2) + 9 * eta2) / 24.0
        let IIIA = nu * F0 * sin(lat) * pow(cos(lat),5) * (61 - 58 * pow(tan(lat),2) + pow(tan(lat),4)) / 720.0
        let IV = nu * F0 * cos(lat)
        let V = nu * F0 * pow(cos(lat),3) * (nu / rho - pow(tan(lat),2)) / 6.0
        let VI = nu * F0 * pow(cos(lat),5) * (5 - 18 * pow(tan(lat),2) + pow(tan(lat),4) + 14 * eta2 - 58 * eta2 * pow(tan(lat),2)) / 120.0

        let northing = I + II * pow((lon - lon0), 2) + III * pow(lon - lon0, 4) + IIIA * pow(lon - lon0, 6)
        let easting = E0 + IV * (lon - lon0) + V * pow(lon - lon0,3) + VI * pow(lon - lon0, 5)

        self.init(easting: easting, northing: northing)
    }
}
