from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 get_moon, Angle, Latitude, Longitude,
                                 UnitSphericalRepresentation)
import datetime
from astropy.extern.six import string_types
import pytz

class Observer(object):

   
    @u.quantity_input(elevation=u.m)
    def __init__(self, location=None, timezone='UTC', name=None, latitude=None,
                 longitude=None, elevation=0*u.m, pressure=None,
                 relative_humidity=None, temperature=None, description=None):
       

        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and
                                 longitude is not None):
            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) an instance of '
                            'astropy.coordinates.EarthLocation or (2) '
                            'latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude.')

        # Accept various timezone inputs, default to UTC
        if isinstance(timezone, datetime.tzinfo):
            self.timezone = timezone
        elif isinstance(timezone, string_types):
            self.timezone = pytz.timezone(timezone)
        else:
            raise TypeError('timezone keyword should be a string, or an '
                            'instance of datetime.tzinfo')

    def __repr__(self):
        
        class_name = self.__class__.__name__
        attr_names = ['name', 'location', 'timezone', 'pressure', 'temperature',
                      'relative_humidity']
        attr_values = [getattr(self, attr) for attr in attr_names]
        attributes_strings = []
        for name, value in zip(attr_names, attr_values):
            if value is not None:
                # Format location for easy readability
                if name == 'location':
                    formatted_loc = ["{} {}".format(i.value, i.unit)
                                     for i in value.to_geodetic()]
                    attributes_strings.append(
                        "{} (lon, lat, el)=({})".format(name,
                                                        ", ".join(formatted_loc)))
                else:
                    if name != 'name':
                        value = repr(value)
                    else:
                        value = "'{}'".format(value)
                    attributes_strings.append("{}={}".format(name, value))
        return "<{}: {}>".format(class_name, ",\n    ".join(attributes_strings))

    @classmethod
    def at_site(cls, site_name, **kwargs):
       
        name = kwargs.pop('name', site_name)
        if 'location' in kwargs:
            raise ValueError("Location kwarg should not be used if "
                             "initializing an Observer with Observer.at_site()")
        return cls(location=EarthLocation.of_site(site_name), name=name, **kwargs)

    def astropy_time_to_datetime(self, astropy_time):
        

        if not astropy_time.isscalar:
            return [self.astropy_time_to_datetime(t) for t in astropy_time]

        # Convert astropy.time.Time to a UTC localized datetime (aware)
        utc_datetime = pytz.utc.localize(astropy_time.utc.datetime)

        # Convert UTC to local timezone
        return self.timezone.normalize(utc_datetime)

    def datetime_to_astropy_time(self, date_time):
        

        if hasattr(date_time, '__iter__'):
            return Time([self.datetime_to_astropy_time(t) for t in date_time])

        # For timezone-naive datetimes, assign local timezone
        if date_time.tzinfo is None:
            date_time = self.timezone.localize(date_time)

        return Time(date_time, location=self.location)

    def _is_broadcastable(self, shp1, shp2):
        """Test if two shape tuples are broadcastable"""
        if shp1 == shp2:
            return True
        for a, b in zip(shp1[::-1], shp2[::-1]):
            if a == 1 or b == 1 or a == b:
                pass
            else:
                return False
        return True

    def _preprocess_inputs(self, time, target=None, grid=True):
        
        if not isinstance(time, Time):
            time = Time(time)

        if target is None:
            return time, None

        # convert any kind of target argument to non-scalar SkyCoord
        target = get_skycoord(target)

        if grid:
            # now we broadcast the targets array so that the first index
            # iterates over targets, any other indices over times
            if not target.isscalar:
                if time.isscalar:
                    target = target[:, np.newaxis]
                while target.ndim <= time.ndim:
                    target = target[:, np.newaxis]
        if not self._is_broadcastable(target.shape, time.shape):
            raise ValueError(
                'Time and Target arguments cannot be broadcast against each other with shapes {} and {}'.format(
                    time.shape, target.shape
                ))
        return time, target

    @profile
    def altaz(self, time, target=None, obswl=None, grid=True):
       
        if target is not None:
            time, target = self._preprocess_inputs(time, target, grid)

        altaz_frame = AltAz(location=self.location, obstime=time,pressure=self.pressure, obswl=obswl, temperature=self.temperature, relative_humidity=self.relative_humidity)
        if target is None:
            # Return just the frame
            return altaz_frame
        else:
            return target.transform_to(altaz_frame)

    def parallactic_angle(self, time, target):
        
        time, coordinate = self._preprocess_inputs(time, target)

        # Eqn (14.1) of Meeus' Astronomical Algorithms
        LST = time.sidereal_time('mean', longitude=self.location.longitude)
        H = (LST - coordinate.ra).radian
        q = np.arctan(np.sin(H) /
                      (np.tan(self.location.latitude.radian) *
                       np.cos(coordinate.dec.radian) -
                       np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad

        return Angle(q)

    # Sun-related methods.
    @u.quantity_input(horizon=u.deg)
    def _horiz_cross(self, t, alt, rise_set, horizon=0*u.degree):
        
        finesse_time_indexes = False
        if alt.ndim == 1:
            raise ValueError('Must supply more at least a 2D grid of altitudes')
        elif alt.ndim == 2:
            # TODO: this test for ndim=2 doesn't work. if times is e.g (2,5)
            # then alt will have ndim=3, but shape (100, 2, 5) so grid
            # is in first index...
            ntargets = alt.shape[1]
            ngrid = alt.shape[0]
            unit = alt.unit
            alt = broadcast_to(alt, (ntargets, ngrid, ntargets)).T
            alt = alt*unit
            extra_dimension_added = True
            if t.shape[1] == 1:
                finesse_time_indexes = True
        else:
            extra_dimension_added = False
        output_shape = (alt.shape[0],) + alt.shape[2:]

        if rise_set == 'rising':
            # Find index where altitude goes from below to above horizon
            condition = (alt[:, :-1, ...] < horizon) * (alt[:, 1:, ...] > horizon)
        elif rise_set == 'setting':
            # Find index where altitude goes from above to below horizon
            condition = (alt[:, :-1, ...] > horizon) * (alt[:, 1:, ...] < horizon)

        noncrossing_indices = np.sum(condition, axis=1, dtype=np.intp) < 1
        alt_lims1 = u.Quantity(np.zeros(output_shape), unit=u.deg)
        alt_lims2 = u.Quantity(np.zeros(output_shape), unit=u.deg)
        jd_lims1 = np.zeros(output_shape)
        jd_lims2 = np.zeros(output_shape)
        if np.any(noncrossing_indices):
            for target_index in set(np.where(noncrossing_indices)[0]):
                warnmsg = ('Target with index {} does not cross horizon={} within '
                           '24 hours'.format(target_index, horizon))
                if (alt[target_index, ...] > horizon).all():
                    warnings.warn(warnmsg, TargetAlwaysUpWarning)
                else:
                    warnings.warn(warnmsg, TargetNeverUpWarning)

            alt_lims1[np.nonzero(noncrossing_indices)] = np.nan
            alt_lims2[np.nonzero(noncrossing_indices)] = np.nan
            jd_lims1[np.nonzero(noncrossing_indices)] = np.nan
            jd_lims2[np.nonzero(noncrossing_indices)] = np.nan

        before_indices = np.array(np.nonzero(condition))
        # we want to add an vector like (0, 1, ...) to get after indices
        array_to_add = np.zeros(before_indices.shape[0])[:, np.newaxis].astype(int)
        array_to_add[1] = 1
        after_indices = before_indices + array_to_add

        al1 = alt[tuple(before_indices)]
        al2 = alt[tuple(after_indices)]
        # slice the time in the same way, but delete the object index
        before_time_index_tuple = np.delete(before_indices, 0, 0)
        after_time_index_tuple = np.delete(after_indices, 0, 0)
        if finesse_time_indexes:
            before_time_index_tuple[1:] = 0
            after_time_index_tuple[1:] = 0
        tl1 = t[tuple(before_time_index_tuple)]
        tl2 = t[tuple(after_time_index_tuple)]

        alt_lims1[tuple(np.delete(before_indices, 1, 0))] = al1
        alt_lims2[tuple(np.delete(before_indices, 1, 0))] = al2
        jd_lims1[tuple(np.delete(before_indices, 1, 0))] = tl1.utc.jd
        jd_lims2[tuple(np.delete(before_indices, 1, 0))] = tl2.utc.jd

        if extra_dimension_added:
            return (alt_lims1.diagonal(), alt_lims2.diagonal(),
                    jd_lims1.diagonal(), jd_lims2.diagonal())
        else:
            return alt_lims1, alt_lims2, jd_lims1, jd_lims2

    @u.quantity_input(horizon=u.deg)
    def _two_point_interp(self, jd_before, jd_after,
                          alt_before, alt_after, horizon=0*u.deg):
        
        slope = (alt_after-alt_before)/((jd_after - jd_before)*u.d)
        crossing_jd = (jd_after*u.d - ((alt_after - horizon)/slope))
        crossing_jd[np.isnan(crossing_jd)] = u.d*MAGIC_TIME.jd
        return np.squeeze(Time(crossing_jd, format='jd'))

    def _altitude_trig(self, LST, target):
        
        LST, target = self._preprocess_inputs(LST, target)
        alt = np.arcsin(np.sin(self.location.latitude.radian) *
                        np.sin(target.dec) +
                        np.cos(self.location.latitude.radian) *
                        np.cos(target.dec) *
                        np.cos(LST.radian - target.ra.radian))
        return alt

    def _calc_riseset(self, time, target, prev_next, rise_set, horizon, N=150):
       
        if not isinstance(time, Time):
            time = Time(time)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altaz = self.altaz(times, target, grid=True)
        altitudes = altaz.alt

        al1, al2, jd1, jd2 = self._horiz_cross(times, altitudes, rise_set,
                                               horizon)
        return self._two_point_interp(jd1, jd2, al1, al2,
                                      horizon=horizon)

    def _calc_transit(self, time, target, prev_next, antitransit=False, N=150):
       
        # TODO FIX BROADCASTING HERE
        if not isinstance(time, Time):
            time = Time(time)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N, for_deriv=True)
        else:
            times = _generate_24hr_grid(time, -1, 0, N, for_deriv=True)

        # The derivative of the altitude with respect to time is increasing
        # from negative to positive values at the anti-transit of the meridian
        if antitransit:
            rise_set = 'rising'
        else:
            rise_set = 'setting'

        altaz = self.altaz(times, target, grid=True)
        altitudes = altaz.alt
        if altitudes.ndim > 2:
            # shape is (M, N, ...) where M is targets and N is grid
            d_altitudes = altitudes.diff(axis=1)
        else:
            # shape is (N, M) where M is targets and N is grid
            d_altitudes = altitudes.diff(axis=0)

        dt = Time((times.jd[1:] + times.jd[:-1])/2, format='jd')

        horizon = 0*u.degree  # Find when derivative passes through zero
        al1, al2, jd1, jd2 = self._horiz_cross(dt, d_altitudes,
                                               rise_set, horizon)
        return self._two_point_interp(jd1, jd2, al1, al2,
                                      horizon=horizon)

    def _determine_which_event(self, function, args_dict):
        
        time = args_dict.pop('time', None)
        target = args_dict.pop('target', None)
        which = args_dict.pop('which', None)
        horizon = args_dict.pop('horizon', None)
        rise_set = args_dict.pop('rise_set', None)
        antitransit = args_dict.pop('antitransit', None)

        # Assemble arguments for function, depending on the function.
        if function == self._calc_riseset:
            args = lambda w: (time, target, w, rise_set, horizon)
        elif function == self._calc_transit:
            args = lambda w: (time, target, w, antitransit)
        else:
            raise ValueError('Function {} not supported in '
                             '_determine_which_event.'.format(function))

        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_event = function(*args('next'))
            if which == 'next':
                return next_event

        if which == 'previous' or which == 'nearest':
            previous_event = function(*args('previous'))
            if which == 'previous':
                return previous_event

        if which == 'nearest':
            mask = abs(time - previous_event) < abs(time - next_event)
            return Time(np.where(mask, previous_event.utc.jd,
                        next_event.utc.jd), format='jd')


        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def target_rise_time(self, time, target, which='nearest', horizon=0*u.degree):
        
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='rising',
                                                horizon=horizon))

    @u.quantity_input(horizon=u.deg)
    def target_set_time(self, time, target, which='nearest', horizon=0*u.degree):
       
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='setting',
                                                horizon=horizon))

    def target_meridian_transit_time(self, time, target, which='nearest'):
        
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which,
                                                rise_set='setting'))

    def target_meridian_antitransit_time(self, time, target, which='nearest'):
       
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which, antitransit=True,
                                                rise_set='setting'))

    @u.quantity_input(horizon=u.deg)
    def sun_rise_time(self, time, which='nearest', horizon=0*u.degree):
       
        return self.target_rise_time(time, get_sun(time), which, horizon)

    @u.quantity_input(horizon=u.deg)
    def sun_set_time(self, time, which='nearest', horizon=0*u.degree):
       
        return self.target_set_time(time, get_sun(time), which, horizon)

    def noon(self, time, which='nearest'):
       
        return self.target_meridian_transit_time(time, get_sun(time), which)

    def midnight(self, time, which='nearest'):
        
        return self.target_meridian_antitransit_time(time, get_sun(time), which)

    # Twilight convenience functions

    def twilight_evening_astronomical(self, time, which='nearest'):
       
        return self.sun_set_time(time, which, horizon=-18*u.degree)

    def twilight_evening_nautical(self, time, which='nearest'):
        
        return self.sun_set_time(time, which, horizon=-12*u.degree)

    def twilight_evening_civil(self, time, which='nearest'):
       
        return self.sun_set_time(time, which, horizon=-6*u.degree)

    def twilight_morning_astronomical(self, time, which='nearest'):
       
        return self.sun_rise_time(time, which, horizon=-18*u.degree)

    def twilight_morning_nautical(self, time, which='nearest'):
       
        return self.sun_rise_time(time, which, horizon=-12*u.degree)

    def twilight_morning_civil(self, time, which='nearest'):
        
        return self.sun_rise_time(time, which, horizon=-6*u.degree)

    # Moon-related methods.

    def moon_rise_time(self, time, **kwargs):
        
        raise NotImplementedError()

    def moon_set_time(self, time, **kwargs):
        
        raise NotImplementedError()

    def moon_illumination(self, time):
        
        if not isinstance(time, Time):
            time = Time(time)

        return moon_illumination(time)

    def moon_phase(self, time=None):
        
        if time is not None and not isinstance(time, Time):
            time = Time(time)

        return moon_phase_angle(time)

    def moon_altaz(self, time, ephemeris=None):
        
        if not isinstance(time, Time):
            time = Time(time)

        moon = get_moon(time, location=self.location, ephemeris=ephemeris)
        return self.altaz(time, moon, grid=False)

    @u.quantity_input(horizon=u.deg)
    def target_is_up(self, time, target, horizon=0*u.degree, return_altaz=False):
       
        if not isinstance(time, Time):
            time = Time(time)

        altaz = self.altaz(time, target)
        observable = altaz.alt > horizon
        if altaz.isscalar:
            observable = bool(observable)
        else:
            # TODO: simply return observable if we move to
            # a fully broadcasted API
            observable = [value for value in observable.flat]

        if not return_altaz:
            return observable
        else:
            return observable, altaz

    @u.quantity_input(horizon=u.deg)
    def is_night(self, time, horizon=0*u.deg, obswl=None):
        
        if not isinstance(time, Time):
            time = Time(time)

        solar_altitude = self.altaz(time, target=get_sun(time), obswl=obswl).alt
        if solar_altitude.isscalar:
            return bool(solar_altitude < horizon)
        else:
            # TODO: simply return solar_altitude < horizon if we move to
            # a fully broadcasted API
            return [val for val in (solar_altitude < horizon).flat]

    def local_sidereal_time(self, time, kind='apparent', model=None):
       
        if not isinstance(time, Time):
            time = Time(time)

        return time.sidereal_time(kind, longitude=self.location.longitude,
                                  model=model)

    def target_hour_angle(self, time, target):
        
        time, target = self._preprocess_inputs(time, target)
        return Longitude(self.local_sidereal_time(time) - target.ra)

    @u.quantity_input(horizon=u.degree)
    def tonight(self, time=None, horizon=0 * u.degree, obswl=None):
        
        current_time = Time.now() if time is None else time
        night_mask = self.is_night(current_time, horizon=horizon, obswl=obswl)
        sun_set_time = self.sun_set_time(current_time, which='next', horizon=horizon)
        # workaround for NPY <= 1.8, otherwise np.where works even in scalar case
        if current_time.isscalar:
            start_time = current_time if night_mask else sun_set_time
        else:
            start_time = np.where(night_mask, current_time, sun_set_time)
            # np.where gives us a list of start Times - convert to Time object
            if not isinstance(start_time, Time):
                start_time = Time(start_time)
        end_time = self.sun_rise_time(start_time, which='next', horizon=horizon)

        return start_time, end_time



apo = Observer.at_site('apo')
time = Time('2001-02-03 04:05:06')
target = SkyCoord(0*u.deg, 0*u.deg)
altaz_frame = apo.altaz(time)
print(altaz_frame)
