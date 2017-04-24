------------------------------------------------------------------------------------------------

-- Max Gold
-- 4/7/2017
-- Script to create multiz46way placental data as postgres database

------------------------------------------------------------------------------------------------

DROP TABLE IF EXISTS multiz46way_placental;

CREATE TABLE IF NOT EXISTS multiz46way_placental(
  ZScore      float DEFAULT NULL,
  FDR_lower   float DEFAULT NULL,
  name        varchar(1024) DEFAULT NULL,
  orientation int DEFAULT NULL,
  chromosome  varchar(1024) DEFAULT NULL,
  LOD         float DEFAULT NULL,
  strand      varchar(1024) DEFAULT NULL,
  "start"     int DEFAULT NULL,
  realhits    float DEFAULT NULL,
  cid         SERIAL PRIMARY KEY,
  FDR         float DEFAULT NULL,
  NLOD        float DEFAULT NULL,
  BBLS        float DEFAULT NULL,
  stop        int DEFAULT NULL,
  medianhits  float DEFAULT NULL,
  accession   varchar(1024) DEFAULT NULL,
  FDR_upper   float DEFAULT NULL,
  BLS         float DEFAULT NULL,
  stdevhits   float DEFAULT NULL
);

\COPY multiz46way_placental FROM '/Users/adaigle/GarNet/data/motifmap/hg19/multiz46way_placental.txt';
