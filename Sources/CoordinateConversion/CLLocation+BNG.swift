import CoreLocation

// A Swift reimplementation of the Objective-C code listed here :
// http://www.hannahfry.co.uk/blog/2014/12/26/more-on-converting-british-national-grid-to-latitude-and-longitude

public extension CLLocation {
    convenience init(_ bng: BNGLocation) {
        self.init(easting: bng.easting, northing: bng.northing)
    }

    convenience init(easting E: Double, northing N: Double) {
        // The Airy 180 semi-major and semi-minor axes used for OSGB36 (m)
        let (a, b) = (6377563.396, 6356256.909)

        let F0 = 0.9996012717        // scale factor on the central meridian
        let lat0 = 49 * Double.pi / 180   // Latitude of true origin (radians)
        let lon0 = -2 * Double.pi / 180   // Longitude of true origin and central meridian (radians)

        // Northing & easting of true origin (m)
        let (N0, E0) = (-100000.0, 400000.0)

        let e2 = 1 - (b * b) / (a * a)     // eccentricity squared
        let n = (a - b) / (a + b)

        // Initialise the iterative variables
        var lat = lat0
        var M = 0.0

        while N - N0 - M >= 0.00001     // Accurate to 0.01mm
        {
            lat = (N - N0 - M) / (a * F0) + lat
            let M1 = (1.0 + n + (5 / 4) * pow(n,2) + (5 / 4) * pow(n,3)) * (lat - lat0)
            let M2 = (3.0 * n + 3.0 * pow(n,2) + (21 / 8) * pow(n,3)) * sin(lat - lat0) * cos(lat + lat0)
            let M3 = ((15 / 8) * pow(n,2) + (15 / 8) * pow(n,3)) * sin(2 * (lat - lat0)) * cos(2 * (lat + lat0))
            let M4 = (35 / 24) * pow(n,3) * sin(3 * (lat - lat0)) * cos(3 * (lat + lat0))

            // meridional arc
            M = b * F0 * (M1 - M2 + M3 - M4)
        }

        // transverse radius of curvature
        let nu = a * F0 / sqrt(1 - e2 * pow(sin(lat),2))

        // meridional radius of curvature
        let rho = a * F0 * (1 - e2) * pow((1 - e2 * pow(sin(lat),2)),(-1.5))
        let eta2 = nu / rho - 1.0

        let secLat = 1.0 / cos(lat)
        let VII = tan(lat) / (2 * rho * nu)

        let VIII = tan(lat) / (24 * rho * pow(nu,3)) * (5 + 3 * pow(tan(lat),2) + eta2 - 9 * pow(tan(lat),2) * eta2)
        let IX = tan(lat) / (720 * rho * pow(nu,5)) * (61 + 90 * pow(tan(lat),2) + 45 * pow(tan(lat),4))
        let X = secLat / nu
        let XI = secLat / (6 * pow(nu,3)) * (nu / rho + 2 * pow(tan(lat),2))
        let XII = secLat / (120 * pow(nu,5)) * (5 + 28 * pow(tan(lat),2) + 24 * pow(tan(lat),4))
        let XIIA = secLat / (5040 * pow(nu,7)) * (61 + 662 * pow(tan(lat),2) + 1320 * pow(tan(lat),4) + 720 * pow(tan(lat),6))
        let dE = E-E0

        // These are on the wrong ellipsoid currently: Airy1830. (Denoted by _1)
        let lat_1 = lat - VII * pow(dE,2) + VIII * pow(dE,4) - IX * pow(dE,6)
        let lon_1 = lon0 + X * dE - XI * pow(dE,3) + XII * pow(dE,5) - XIIA * pow(dE,7)

        // Want to convert to the GRS80 ellipsoid.
        // First convert to cartesian from spherical polar coordinates
        var H = 0.0     // Third spherical coord.
        let x_1 = (nu / F0 + H) * cos(lat_1) * cos(lon_1)
        let y_1 = (nu / F0 + H) * cos(lat_1) * sin(lon_1)
        let z_1 = ((1 - e2) * nu / F0 + H) * sin(lat_1)

        // Perform Helmut transform (to go between Airy 1830 (_1) and GRS80 (_2))
        let s = -20.4894 * pow(10.0,-6)           // The scale factor -1
        let tx = 446.448                        // The translations along x,y,z axes respectively
        let ty = -125.157
        let tz = 542.060
        let rxs = 0.1502                        // The rotations along x,y,z respectively, in seconds
        let rys = 0.2470
        let rzs = 0.8421
        let rx = rxs * Double.pi / (180 * 3600.0)    // In radians
        let ry = rys * Double.pi / (180 * 3600.0)
        let rz = rzs * Double.pi / (180 * 3600.0)
        let x_2 = tx + (1 + s) * x_1 + (-rz) * y_1 + (ry) * z_1
        let y_2 = ty + (rz) * x_1 + (1 + s) * y_1 + (-rx) * z_1
        let z_2 = tz + (-ry) * x_1 + (rx) * y_1 + (1 + s) * z_1

        // Back to spherical polar coordinates from cartesian
        // Need some of the characteristics of the new ellipsoid
        let a_2 = 6378137.000                   // The GSR80 semi-major and semi-minor axes used for WGS84(m)
        let b_2 = 6356752.3141
        let e2_2 = 1 - (b_2 * b_2) / (a_2 * a_2) // The eccentricity of the GRS80 ellipsoid
        let p = sqrt(pow(x_2,2) + pow(y_2,2))

        // Lat is obtained by an iterative proceedure:
        lat = atan2(z_2,(p * (1 - e2_2)))       // Initial value
        var latold = 2 * Double.pi

        var nu_2 = 0.0

        while abs(lat - latold) > pow(10,-16)
        {
            let latTemp = lat
            lat = latold
            latold = latTemp
            nu_2 = a_2 / sqrt(1 - e2_2 * pow(sin(latold),2))
            lat = atan2(z_2 + e2_2 * nu_2 * sin(latold), p)
        }

        // Lon and height are then pretty easy
        var lon = atan2(y_2, x_2)
        H = p / cos(lat) - nu_2

        // Convert to degrees
        lat = lat * 180 / Double.pi
        lon = lon * 180 / Double.pi

        self.init(latitude: lat, longitude: lon)
    }

    func toBNG() -> BNGLocation {
        BNGLocation(self)
    }
}
