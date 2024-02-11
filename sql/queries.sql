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

/* Counts the number of planets in each solar system where
the corresponding stars are larger than the Sun. Return the star radius
and the number of planets in the system only where n > 1. Sort
in descending order based on star radii */
select s.radius, count(p.koi_name)
from Star as s 
join Planet as p using (kepler_id)
where s.radius > 1
group by s.kepler_id
having count(p.koi_name) > 1
order by s.radius desc

/* Return ID, temperature, and radius for all stars which do not
have a planetary companion. Sort by temperature in descending order. */
select s.kepler_id, s.t_eff, s.radius
from Star as s
full outer join Planet as p using (kepler_id)
where p.koi_name is null
order by s.t_eff desc 

/* Return average planet equilibrium temperature (rounded to 1 decimal place), 
and min and max star effective temperature. */
select round(avg(p.t_eq), 1), min(s.t_eff), max(s.t_eff)
from Star as s
join Planet as p using (kepler_id)
where s.t_eff > (
  select avg(t_eff) from Star
  );

/* Return koi_name and radius of planets with orbit the 5 largest stars. 
Return the radius of the corresponding star. */
select p.koi_name, p.radius, s.radius
from Star as s
join Planet as p using (kepler_id)
where s.kepler_id in (
  select kepler_id
  from Star
  order by radius desc
  limit 5
);

/* Update Kepler names of planets without a confirmed status to NULL.
Delete rows with negative radius. */
select status from Planet;

update Planet
set kepler_name = NULL
where status = 'FALSE POSITIVE' or status = 'CANDIDATE';

delete from Planet
where radius < 0;

/* Set up Planet table and fill it with 3 rows. */
create table Planet (
  kepler_id integer NOT NULL,
  koi_name varchar(15) NOT NULL UNIQUE,
  kepler_name varchar(15),
  status varchar(20) NOT NULL,
  radius float NOT NULL
);

insert into Planet (kepler_id, koi_name, kepler_name, status, radius)
values (6862328, 'K00865.01', NULL, 'CANDIDATE', 119.021), 
       (10187017, 'K00082.05', 'Kepler-102 b', 'CONFIRMED', 5.286),
       (10187017, 'K00082.04', 'Kepler-102 c', 'CONFIRMED', 7.071);

select * from Planet;

/* Create the Star and Planet tables and fill with data from existing 
CSV files. Link Star and Planet kepler_id. */
create table Star (
  kepler_id integer PRIMARY KEY,
  t_eff integer NOT NULL,
  radius float NOT NULL
 );
 
copy Star (kepler_id, t_eff, radius) from 'stars.csv' csv; 

create table Planet (
  kepler_id integer REFERENCES Star (kepler_id),
  koi_name varchar(20) PRIMARY KEY,
  kepler_name varchar(20),
  status varchar(20) NOT NULL,
  period float,
  radius float,
  t_eq integer
);

copy Planet (kepler_id, koi_name, kepler_name, status, period, radius, t_eq)
from 'planets.csv' csv;

/* Add 2 columns to existing Star table for equatorial coordinates, then
fill new columns with data. */
alter table Star
add column ra FLOAT,
add column decl FLOAT;

delete from Star;

copy Star (kepler_id, t_eff, radius, ra, decl) from 'stars_full.csv' csv;
select * from Star;


