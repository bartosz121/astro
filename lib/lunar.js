"use strict";
const julian = require("./julian");
const sidereal = require("./sidereal");

const RAD = Math.PI / 180;
const AU = 1 / 149597870.7; // convert km to AU

function limit_angle(radians) {
  return radians - Math.floor(radians / (2 * Math.PI)) * (2 * Math.PI);
}

/* Lunar position, accurate to about ten minutes between 1950-2050,
 * http://aa.quae.nl/en/reken/hemelpositie.html#2 */
const Q0 = 218.316 * RAD;
const Q1 = 13.176396 * RAD;
const M0 = 134.963 * RAD;
const M1 = 13.064993 * RAD;
const F0 = 93.272 * RAD;
const F1 = 13.22935 * RAD;
const L0 = 6.289 * RAD;
const L1 = 5.128 * RAD;
const E0 = 23.439 * RAD;
const E1 = -0.00000036 * RAD;
const D0 = 385001 * AU;
const D1 = -20905 * AU;

/* Geocentric ecliptic longitude compared to the equinox */
function _ecliptic_equinox_longitude(t) {
  return Q0 + Q1 * t;
}

function _mean_anomaly(t) {
  return M0 + M1 * t;
}

/* Mean distance of the moon to its ascending node */
function _mean_distance(t) {
  return F0 + F1 * t;
}

/* Ecliptic longitude */
function _longitude(t) {
  return _ecliptic_equinox_longitude(t) + L0 * Math.sin(_mean_anomaly(t));
}

/* Ecliptic latitude */
function _latitude(t) {
  return L1 * Math.sin(_mean_distance(t));
}

function _obliquity(t) {
  return E0 + E1 * t;
}

function _distance(t) {
  return D0 + D1 * Math.cos(_mean_anomaly(t));
}

function _right_ascension_direct(lat, lon, obl) {
  return Math.atan2(
    Math.sin(lon) * Math.cos(obl) - Math.tan(lat) * Math.sin(obl),
    Math.cos(lon)
  );
}

function _right_ascension(t) {
  return _right_ascension_direct(_latitude(t), _longitude(t), _obliquity(t));
}

function _declination_direct(lat, lon, obl) {
  return Math.asin(
    Math.sin(lat) * Math.cos(obl) +
      Math.cos(lat) * Math.sin(obl) * Math.sin(lon)
  );
}

function _declination(t) {
  return _declination_direct(_latitude(t), _longitude(t), _obliquity(t));
}

function _hour_angle_direct(lst, ra) {
  return limit_angle(lst) - limit_angle(ra);
}

function _transit_direct(t, ha) {
  return t - ha * (12 / Math.PI / 24);
}

function _hour_angle(t, lon) {
  return _hour_angle_direct(sidereal.local(t, lon), _right_ascension(t));
}

function _hour_angle_iterative(t, lon, ha0) {
  // Iteratively improve by walking back in time toward a smaller hour angle.
  const t1 = _transit_direct(t, ha0);
  const ha1 = _hour_angle(t1, lon) + ha0;
  return ha1;
}

function _hour_angle_refined(t, lon) {
  let ha = _hour_angle(t, lon);

  ha = _hour_angle_iterative(t, lon, ha);
  ha = _hour_angle_iterative(t, lon, ha);
  ha = _hour_angle_iterative(t, lon, ha);
  ha = _hour_angle_iterative(t, lon, ha);

  return ha;
}

function _hours_later(t, hrs) {
  return t + hrs / 24;
}

function _transit(t, lat, lon) {
  // Go 24h ahead and look backwards for the transit...
  const later_t = _hours_later(t, 24);
  let result = _transit_direct(later_t, _hour_angle_refined(later_t, lon));
  if (result <= later_t) {
    return result;
  }
  // If we've passed our max time though... just look near the beginning
  result = _transit_direct(t, _hour_angle_refined(t, lon));
  if (result >= t) {
    return result;
  }
  return NaN;
}

function _elevation(t, lat, lon) {
  const decl = _declination(t);
  const ha = _hour_angle(t, lon);

  return Math.asin(
    Math.sin(lat) * Math.sin(decl) +
      Math.cos(lat) * Math.cos(decl) * Math.cos(ha)
  );
}

/* http://www.stargazing.net/kepler/moonrise.html article */
const PARALLAX = 0.0023212879051524586;

function _rise_and_set(t, lat, lon) {
  const h = -PARALLAX;
  let h0 = _elevation(t, lat, lon) - h;
  let h1 = h0;
  // Start at the beginning of the day
  let moonrise = NaN;
  let moonset = NaN;
  // Go in 2 hour chunks.
  for (let i = 0; i <= 24; i += 2) {
    if (i !== 0) {
      h1 = _elevation(_hours_later(t, i), lat, lon) - h;
    }
    const h2 = _elevation(_hours_later(t, i + 1), lat, lon) - h;

    // Fit h0, h1, h2 to a parabola
    const a = (h2 + h0) / 2 - h1;
    const b = (h2 - h0) / 2;
    const xe = -b / (2 * a); // vertex of parabola
    const ye = (a * xe + b) * xe + h1;

    // Discriminant
    const d = b * b - 4 * a * h1;
    let roots = 0;
    let x1;
    let x2;

    // Count roots
    if (d >= 0) {
      const dx = Math.sqrt(d) / (Math.abs(a) * 2);
      x1 = xe - dx;
      x2 = xe + dx;

      if (Math.abs(x1) <= 1) {
        roots++;
      }
      if (Math.abs(x2) <= 1) {
        roots++;
      }
      if (x1 < -1) {
        x1 = x2;
      }
    }

    if (roots === 1) {
      if (h0 < 0 && isNaN(moonrise)) {
        moonrise = i + x1;
      } else if (isNaN(moonset)) {
        moonset = i + x1;
      }
    } else if (roots === 2) {
      if (isNaN(moonrise)) {
        moonrise = i + (ye < 0 ? x2 : x1);
      }
      if (isNaN(moonset)) {
        moonset = i + (ye < 0 ? x1 : x2);
      }
    }

    // Found all the things!
    if (moonrise < 24 && moonset < 24) {
      break;
    }
    // Move two hours of elevation
    h0 = h2;
  }
  return {
    moonrise: _hours_later(t, moonrise),
    moonset: _hours_later(t, moonset),
  };
}

// http://www.ben-daglish.net/moon.shtml
function _phaseNumber(ts) {
  const date = _dateFromTimestamp(ts);
  const year = date.getFullYear();
  const month = date.getMonth() + 1; // January == 0
  const day = date.getDate();
  let r = year % 100;
  r %= 19;
  if (r > 9) {
    r -= 19;
  }
  r = ((r * 11) % 30) + parseInt(month) + parseInt(day);
  if (month < 3) {
    r += 2;
  }
  r -= year < 2000 ? 4 : 8.3;
  r = Math.floor(r + 0.5) % 30;
  return r < 0 ? r + 30 : r;
  // http://www.stargazing.net/david/moon/index29days.html
}

function _dayFromTimestamp(ts) {
  return ts / 86400000 + 2440587.5;
}

function _dateFromTimestamp(ts) {
  if (typeof ts === "string") {
    ts = parseInt(ts);
  }
  return new Date(ts);
}

function _phaseString(phaseNumber) {
  let moonPhase = null;
  phaseNumber = Math.round(phaseNumber);
  console.log(phaseNumber);

  // New Moon = 0
  // Waxing Crescent = 1-7
  // First Quarter = 8
  // Waxing Gibbous = 9-14
  // Full Moon = 15
  // Waning Gibbous = 16-21
  // Last Quarter = 22
  // Waning Crescent = 22-29
  switch (true) {
    case phaseNumber === 0:
      moonPhase = "new-moon";
      break;
    case 0 < phaseNumber && phaseNumber < 8:
      moonPhase = "Waxing-crescent";
      break;
    case phaseNumber === 8:
      moonPhase = "first-quarter";
      break;
    case 8 < phaseNumber && phaseNumber < 15:
      moonPhase = "waxing-gibbous";
      break;
    case phaseNumber === 15:
      moonPhase = "full-moon";
      break;
    case 15 < phaseNumber && phaseNumber < 22:
      moonPhase = "waning-gibbous";
      break;
    case phaseNumber === 22:
      moonPhase = "last-quarter";
      break;
    case 22 < phaseNumber && phaseNumber < 31:
      moonPhase = "waning-crescent";
      break;
  }

  return moonPhase;
}

function latitude(ms) {
  return _latitude(julian.to(ms));
}

function longitude(ms) {
  return _longitude(julian.to(ms));
}

function obliquity(ms) {
  return _obliquity(julian.to(ms));
}

function distance(ms) {
  return _distance(julian.to(ms));
}

function declination(ms) {
  return _declination(julian.to(ms));
}

function right_ascension(ms) {
  return _right_ascension(julian.to(ms));
}

function elevation(ms, lat, lon) {
  return _elevation(julian.to(ms), lat * RAD, lon * RAD);
}

function transit(ms, lat, lon) {
  return julian.from(_transit(julian.to(ms), lat * RAD, lon * RAD));
}

// In the next 24 hours
function rise(ms, lat, lon) {
  return julian.from(
    _rise_and_set(julian.to(ms), lat * RAD, lon * RAD).moonrise
  );
}

// In the next 24 hours
function set(ms, lat, lon) {
  return julian.from(
    _rise_and_set(julian.to(ms), lat * RAD, lon * RAD).moonset
  );
}

// Sometimes, you just want everything
function rise_and_set(ms, lat, lon) {
  const x = _rise_and_set(julian.to(ms), lat * RAD, lon * RAD);
  x.moonrise = julian.from(x.moonrise);
  x.moonset = julian.from(x.moonset);
  return x;
}

function phase(ms) {
  const phaseNumber = _phaseNumber(ms);
  const p = {
    number: phaseNumber,
    string: _phaseString(phaseNumber),
  };

  return p;
}

exports.latitude = latitude;
exports.longitude = longitude;
exports.obliquity = obliquity;
exports.distance = distance;
exports.declination = declination;
exports.right_ascension = right_ascension;
exports.elevation = elevation;
exports.transit = transit;
exports.rise = rise;
exports.set = set;
exports.rise_and_set = rise_and_set;
exports.phase = phase;
