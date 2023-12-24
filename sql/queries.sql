/* Check confirmed exoplanets. */
select kepler_name, radius from planet
where kepler_name is not NULL
and radius between 1 and 3

/* Analyse size of unconfirmed exoplanets. */
select MIN(radius), MAX(radius), AVG(radius), STDDEV(radius)
from planet where kepler_name is NULL

/* Find out how many planets are in a multi-planet system. */
select kepler_id, count(koi_name)
from planet
group by kepler_id
having count(kepler_id) > 1
order by count(kepler_id) desc

/* Return radius of star and planet pairs whose radii have a 
ratio greater than the Sun-to-Earth radius ratio. */
select s.radius as sun_radius, p.radius as planet_radius
from Star as s, Planet as p
where s.kepler_id = p.kepler_id and s.radius > p.radius
order by sun_radius desc 