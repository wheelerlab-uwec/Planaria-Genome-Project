Submit-block ::= {
  contact {
    contact {
      name name {
        last "Wheeler",
        first "Nicolas",
        middle "",
        initials "J.",
        suffix "",
        title ""
      },
      affil std {
        affil "University of Wisconsin-Eau Claire",
        div "Department of Biology",
        city "Eau Claire",
        sub "WI",
        country "USA",
        street "105 Garfield Ave",
        email "wheelenj@uwec.edu",
        postal-code "54701"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Wheeler",
            first "Nicolas",
            initials "N.J.",
            suffix ""
          }
        }
      },
      affil std {
        affil "University of Wisconsin-Eau Claire",
        div "Department of Biology",
        city "Eau Claire",
        sub "WI",
        country "USA",
        street "105 Garfield Ave",
        postal-code "54701"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Wheeler",
              first "Nicolas",
              initials "N.J.",
              suffix ""
            }
          }
        }
      },
      title "Genome assembly of the freshwater planarian Girardia dorotocephala"
    }
  }
}
Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA1545025"
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        "SAMN64250158"
      }
    }
  }
}
