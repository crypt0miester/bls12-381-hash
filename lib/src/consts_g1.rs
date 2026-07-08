// BLS12-381 G1 SSWU + iso-11 constants, Montgomery form (R = 2^384 mod p),
// limbs little-endian; extracted from bls12_381 0.8.0 map_g1.rs / fp.rs.

pub const MODULUS: [u64; 6] = [
    0xb9feffffffffaaab,
    0x1eabfffeb153ffff,
    0x6730d2a0f6b0f624,
    0x64774b84f38512bf,
    0x4b1ba7b6434bacd7,
    0x1a0111ea397fe69a,
];

pub const INV: u64 = 0x89f3fffcfffcfffd;

pub const R2: [u64; 6] = [
    0xf4df1f341c341746,
    0x0a76e6a609d104f1,
    0x8de5476c4c95b6d5,
    0x67eb88a9939d83c0,
    0x9a793e85b519952d,
    0x11988fe592cae3aa,
];

pub const R: [u64; 6] = [
    0x760900000002fffd,
    0xebf4000bc40c0002,
    0x5f48985753c758ba,
    0x77ce585370525745,
    0x5c071a97a256ec6d,
    0x15f65ec3fa80e493,
];

pub const SSWU_ELLP_A: [u64; 6] = [
    0x2f65aa0e9af5aa51,
    0x86464c2d1e8416c3,
    0xb85ce591b7bd31e2,
    0x27e11c91b5f24e7c,
    0x28376eda6bfc1835,
    0x155455c3e5071d85,
];

pub const SSWU_ELLP_B: [u64; 6] = [
    0xfb996971fe22a1e0,
    0x9aa93eb35b742d6f,
    0x8c476013de99c5c4,
    0x873e27c3a221e571,
    0xca72b5e45a52d888,
    0x06824061418a386b,
];

pub const SSWU_XI: [u64; 6] = [
    0x886c00000023ffdc,
    0x0f70008d3090001d,
    0x77672417ed5828c3,
    0x9dac23e943dc1740,
    0x50553f1b9c131521,
    0x078c712fbe0ab6e8,
];

pub const ISO11_XNUM: [[u64; 6]; 12] = [
    [0x4d18b6f3af00131c, 0x19fa219793fee28c, 0x3f2885f1467f19ae, 0x23dcea34f2ffb304, 0xd15b58d2ffc00054, 0x0913be200a20bef4],
    [0x898985385cdbbd8b, 0x3c79e43cc7d966aa, 0x1597e193f4cd233a, 0x8637ef1e4d6623ad, 0x11b22deed20d827b, 0x07097bc5998784ad],
    [0xa542583a480b664b, 0xfc7169c026e568c6, 0x5ba2ef314ed8b5a6, 0x5b5491c05102f0e7, 0xdf6e99707d2a0079, 0x0784151ed7605524],
    [0x494e212870f72741, 0xab9be52fbda43021, 0x26f5577994e34c3d, 0x049dfee82aefbd60, 0x65dadd7828505289, 0x0e93d431ea011aeb],
    [0x90ee774bd6a74d45, 0x7ada1c8a41bfb185, 0x0f1a8953b325f464, 0x104c24211be4805c, 0x169139d319ea7a8f, 0x09f20ead8e532bf6],
    [0x6ddd93e2f43626b7, 0xa5482c9aa1ccd7bd, 0x143245631883f4bd, 0x2e0a94ccf77ec0db, 0xb0282d480e56489f, 0x18f4bfcbb4368929],
    [0x23c5f0c953402dfd, 0x7a43ff6958ce4fe9, 0x2c390d3d2da5df63, 0xd0df5c98e1f9d70f, 0xffd89869a572b297, 0x1277ffc72f25e8fe],
    [0x79f4f0490f06a8a6, 0x85f894a88030fd81, 0x12da3054b18b6410, 0xe2a57f6505880d65, 0xbba074f260e400f1, 0x08b76279f621d028],
    [0xe67245ba78d5b00b, 0x8456ba9a1f186475, 0x7888bff6e6b33bb4, 0xe21585b9a30f86cb, 0x05a69cdcef55feee, 0x09e699dd9adfa5ac],
    [0x0de5c357bff57107, 0x0a0db4ae6b1a10b2, 0xe256bb67b3b3cd8d, 0x8ad456574e9db24f, 0x0443915f50fd4179, 0x098c4bf7de8b6375],
    [0xe6b0617e7dd929c7, 0xfe6e37d442537375, 0x1dafdeda137a489e, 0xe4efd1ad3f767ceb, 0x4a51d8667f0fe1cf, 0x054fdf4bbf1d821c],
    [0x72db2a50658d767b, 0x8abf91faa257b3d5, 0xe969d6833764ab47, 0x464170142a1009eb, 0xb14f01aadb30be2f, 0x18ae6a856f40715d],
];

pub const ISO11_XDEN: [[u64; 6]; 11] = [
    [0xb962a077fdb0f945, 0xa6a9740fefda13a0, 0xc14d568c3ed6c544, 0xb43fc37b908b133e, 0x9c0b3ac929599016, 0x0165aa6c93ad115f],
    [0x23279a3ba506c1d9, 0x92cfca0a9465176a, 0x3b294ab13755f0ff, 0x116dda1c5070ae93, 0xed4530924cec2045, 0x083383d6ed81f1ce],
    [0x9885c2a6449fecfc, 0x4a2b54ccd37733f0, 0x17da9ffd8738c142, 0xa0fba72732b3fafd, 0xff364f36e54b6812, 0x0f29c13c660523e2],
    [0xe349cc118278f041, 0xd487228f2f3204fb, 0xc9d325849ade5150, 0x43a92bd69c15c2df, 0x1c2c7844bc417be4, 0x12025184f407440c],
    [0x587f65ae6acb057b, 0x1444ef325140201f, 0xfbf995e71270da49, 0xccda066072436a42, 0x7408904f0f186bb2, 0x13b93c63edf6c015],
    [0xfb918622cd141920, 0x4a4c64423ecaddb4, 0x0beb232927f7fb26, 0x30f94df6f83a3dc2, 0xaeedd424d780f388, 0x06cc402dd594bbeb],
    [0xd41f761151b23f8f, 0x32a92465435719b3, 0x64f436e888c62cb9, 0xdf70a9a1f757c6e4, 0x6933a38d5b594c81, 0x0c6f7f7237b46606],
    [0x693c08747876c8f7, 0x22c9850bf9cf80f0, 0x8e9071dab950c124, 0x89bc62d61c7baf23, 0xbc6be2d8dad57c23, 0x17916987aa14a122],
    [0x1be3ff439c1316fd, 0x9965243a7571dfa7, 0xc7f7f62962f5cd81, 0x32c6aa9af394361c, 0xbbc2ee18e1c227f4, 0x0c102cbac531bb34],
    [0x997614c97bacbf07, 0x61f86372b99192c0, 0x5b8c95fc14353fc3, 0xca2b066c2a87492f, 0x16178f5bbf698711, 0x12a6dcd7f0f4e0e8],
    [0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba, 0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493],
];

pub const ISO11_YNUM: [[u64; 6]; 16] = [
    [0x2b567ff3e2837267, 0x1d4d9e57b958a767, 0xce028fea04bd7373, 0xcc31a30a0b6cd3df, 0x7d7b18a682692693, 0x0d300744d42a0310],
    [0x99c2555fa542493f, 0xfe7f53cc4874f878, 0x5df0608b8f97608a, 0x14e03832052b49c8, 0x706326a6957dd5a4, 0x0a8dadd9c2414555],
    [0x13d942922a5cf63a, 0x357e33e36e261e7d, 0xcf05a27c8456088d, 0x0000bd1de7ba50f0, 0x83d0c7532f8c1fde, 0x13f70bf38bbf2905],
    [0x5c57fd95bfafbdbb, 0x28a359a65e541707, 0x3983ceb4f6360b6d, 0xafe19ff6f97e6d53, 0xb3468f4550192bf7, 0x0bb6cde49d8ba257],
    [0x590b62c7ff8a513f, 0x314b4ce372cacefd, 0x6bef32ce94b8a800, 0x6ddf84a095713d5f, 0x64eace4cb0982191, 0x0386213c651b888d],
    [0xa5310a31111bbcdd, 0xa14ac0f5da148982, 0xf9ad9cc95423d2e9, 0xaa6ec095283ee4a7, 0xcf5b1f022e1c9107, 0x01fddf5aed881793],
    [0x65a572b0d7a7d950, 0xe25c2d8183473a19, 0xc2fcebe7cb877dbd, 0x05b2d36c769a89b0, 0xba12961be86e9efb, 0x07eb1b29c1dfde1f],
    [0x93e09572f7c4cd24, 0x364e929076795091, 0x8569467e68af51b5, 0xa47da89439f5340f, 0xf4fa918082e44d64, 0x0ad52ba3e6695a79],
    [0x911429844e0d5f54, 0xd03f51a3516bb233, 0x3d587e5640536e66, 0xfa86d2a3a9a73482, 0xa90ed5adf1ed5537, 0x149c9c326a5e7393],
    [0x462bbeb03c12921a, 0xdc9af5fa0a274a17, 0x9a558ebde836ebed, 0x649ef8f11a4fae46, 0x8100e1652b3cdc62, 0x1862bd62c291dacb],
    [0x05c9b8ca89f12c26, 0x0194160fa9b9ac4f, 0x6a643d5a6879fa2c, 0x14665bdd8846e19d, 0xbb1d0d53af3ff6bf, 0x12c7e1c3b28962e5],
    [0xb55ebf900b8a3e17, 0xfedc77ec1a9201c4, 0x1f07db10ea1a4df4, 0x0dfbd15dc41a594d, 0x389547f2334a5391, 0x02419f98165871a4],
    [0xb416af000745fc20, 0x8e563e9d1ea6d0f5, 0x7c763e17763a0652, 0x01458ef0159ebbef, 0x8346fe421f96bb13, 0x0d2d7b829ce324d2],
    [0x93096bb538d64615, 0x6f2a2619951d823a, 0x8f66b3ea59514fa4, 0xf563e63704f7092f, 0x724b136c4cf2d9fa, 0x046959cfcfd0bf49],
    [0xea748d4b6e405346, 0x91e9079c2c02d58f, 0x41064965946d9b59, 0xa06731f1d2bbe1ee, 0x07f897e267a33f1b, 0x1017290919210e5f],
    [0x872aa6c17d985097, 0xeecc53161264562a, 0x07afe37afff55002, 0x54759078e5be6838, 0xc4b92d15db8acca8, 0x106d87d1b51d13b9],
];

pub const ISO11_YDEN: [[u64; 6]; 16] = [
    [0xeb6c359d47e52b1c, 0x18ef5f8a10634d60, 0xddfa71a0889d5b7e, 0x723e71dcc5fc1323, 0x52f45700b70d5c69, 0x0a8b981ee47691f1],
    [0x616a3c4f5535b9fb, 0x6f5f037395dbd911, 0xf25f4cc5e35c65da, 0x3e50dffea3c62658, 0x6a33dca523560776, 0x0fadeff77b6bfe3e],
    [0x2be9b66df470059c, 0x24a2c159a3d36742, 0x115dbe7ad10c2a37, 0xb6634a652ee5884d, 0x04fe8bb2b8d81af4, 0x01c2a7a256fe9c41],
    [0xf27bf8ef3b75a386, 0x898b367476c9073f, 0x24482e6b8c2f4e5f, 0xc8e0bbd6fe110806, 0x59b0c17f7631448a, 0x11037cd58b3dbfbd],
    [0x31c7912ea267eec6, 0x1dbf6f1c5fcdb700, 0xd30d4fe3ba86fdb1, 0x3cae528fbee9a2a4, 0xb1cce69b6aa9ad9a, 0x044393bb632d94fb],
    [0xc66ef6efeeb5c7e8, 0x9824c289dd72bb55, 0x71b1a4d2f119981d, 0x104fc1aafb0919cc, 0x0e49df01d942a628, 0x096c3a09773272d4],
    [0x9abc11eb5fadeff4, 0x32dca50a885728f0, 0xfb1fa3721569734c, 0xc4b76271ea6506b3, 0xd466a75599ce728e, 0x0c81d4645f4cb6ed],
    [0x4199f10e5b8be45b, 0xda64e495b1e87930, 0xcb353efe9b33e4ff, 0x9e9efb24aa6424c6, 0xf08d33680a237465, 0x0d3378023e4c7406],
    [0x7eb4ae92ec74d3a5, 0xc341b4aa9fac3497, 0x5be603899e907687, 0x03bfd9cca75cbdeb, 0x564c2935a96bfa93, 0x0ef3c33371e2fdb5],
    [0x7ee91fd449f6ac2e, 0xe5d5bd5cb9357a30, 0x773a8ca5196b1380, 0xd0fda172174ed023, 0x6cb95e0fa776aead, 0x0d22d5a40cec7cff],
    [0xf727e09285fd8519, 0xdc9d55a83017897b, 0x7549d8bd057894ae, 0x178419613d90d8f8, 0xfce95ebdeb5b490a, 0x0467ffaef23fc49e],
    [0xc1769e6a7c385f1b, 0x79bc930deac01c03, 0x5461c75a23ede3b5, 0x6e20829e5c230c45, 0x828e0f1e772a53cd, 0x116aefa749127bff],
    [0x101c10bf2744c10a, 0xbbf18d053a6a3154, 0xa0ecf39ef026f602, 0xfc009d4996dc5153, 0xb9000209d5bd08d3, 0x189e5fe4470cd73c],
    [0x7ebd546ca1575ed2, 0xe47d5a981d081b55, 0x57b2b625b6d4ca21, 0xb0a1ba04228520cc, 0x98738983c2107ff3, 0x13dddbc4799d81d6],
    [0x09319f2e39834935, 0x039e952cbdb05c21, 0x55ba77a9a2f76493, 0xfd04e3dfc6086467, 0xfb95832e7d78742e, 0x0ef9c24eccaf5e0e],
    [0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba, 0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493],
];


pub const SSWU_C1_NEG_B_OVER_A: [u64; 6] = [
    0x052583c93555a7fe,
    0x3b40d72430f93c82,
    0x1b75faa0105ec983,
    0x2527e7dc63851767,
    0x99fffd1f34fc181d,
    0x097cab54770ca0d3,
];

// Shallue-van de Woestijne direct-map constants (Wahby-Boneh eprint
// 2019/403 section 3, u0 = -3), Montgomery form, derived offline.

pub const SVDW_U0: [u64; 6] = [
    0xcbe1fffffff6000a,
    0x9827ffd8c7d7fff7,
    0x17b8aedce8bcd83b,
    0xc5fad9948998326e,
    0xcd3da75be2de413d,
    0x0c201972bcfd0614,
];

pub const SVDW_NEG_U0: [u64; 6] = [
    0xee1d00000009aaa1,
    0x86840025e97c0007,
    0x4f7823c40df41de8,
    0x9e7c71f069ece051,
    0x7dde005a606d6b99,
    0x0de0f8777c82e085,
];

pub const SVDW_F_U0: [u64; 6] = [
    0xed1cffffffb455a1,
    0x3283fed73d7bffc1,
    0x804ac4babeea4207,
    0x15c7f6e3eeff9fb8,
    0x9985b69dac1a42fe,
    0x0ef2e2b0fc697ad0,
];

pub const SVDW_SQRT_M27: [u64; 6] = [
    0x6039bb5b26b74d45,
    0x007ecdd28055e2fa,
    0xf6f4c41d3aecfe90,
    0x9bf50546bf2f358a,
    0xd31102abe389437c,
    0x077b4f5a1dd8a46d,
];

pub const SVDW_C1: [u64; 6] = [
    0x272b5dad93607bf3,
    0x438166fc34e8f181,
    0x233673f0a4708e3c,
    0x1d38bb9b948e0aee,
    0xa877818321fb578b,
    0x0aae23e8cd2dc279,
];

pub const SVDW_INV_3U0SQ: [u64; 6] = [
    0x158e425ed097b74f,
    0x5dadc71c7e2c4bda,
    0x9d5d01ae2fc08e96,
    0x482181f1982a7a90,
    0x2324e6d352d74573,
    0x0884b37c10d55646,
];

pub const SVDW_B: [u64; 6] = [
    0xaa270000000cfff3,
    0x53cc0032fc34000a,
    0x478fe97a6b0a807f,
    0xb1d37ebee6ba24d7,
    0x8ec9733bbf78ab2f,
    0x09d645513d83de7e,
];

// Knuth-adapted evaluation constants for the iso-11 polynomials
// (TAOCP 4.6.4 style preprocessing): each polynomial becomes a chain
// of quadratic stages (w + beta*x + gamma) with the betas baked into
// the generated evaluator as additions. 

// Derived and expansion-checked offline, layout per polynomial: 
// base constants, then one (gamma, eps) pair per stage, 
// then the leading coefficient if not monic.

pub const ISO11A_XNUM: [[u64; 6]; 12] = [
    [
        0xef4f14c97b9fbf14,
        0x0e2c633fbd5d92b5,
        0x13fcac81a92abf44,
        0x185787ad43cd2458,
        0x874e1c1ffff0dbe2,
        0x08d09786b3710269,
    ],
    [
        0xb168cef405b4ac5b,
        0x53740513cb321aa0,
        0x98deb96b0403f062,
        0x9e71b55c9453d135,
        0x19badfac89987403,
        0x0769f22b0ae96e39,
    ],
    [
        0x3e1b9d5863ff2606,
        0x4cad7793ce42ecf6,
        0x39c8bb4e8924ee2d,
        0x0184510d06a499dd,
        0x378c8b62b027d8e4,
        0x0bc98e64778a1c02,
    ],
    [
        0x8826be632f6f3dd7,
        0xc62246e2162da0da,
        0xd5b076e69b1aad5b,
        0x304e66dd8a3ac877,
        0x04b4d25d99be19a6,
        0x104556b8fb8d7775,
    ],
    [
        0xd51465e2886cb022,
        0xea71085f57b0bcf0,
        0x221f813eaf389916,
        0x33192e18f4a1f0d3,
        0xfd2b076d0399d555,
        0x090f8a95de53631f,
    ],
    [
        0xb4ddd1107d398a5b,
        0x096a86fddd005e03,
        0xe0bfe0d6132c71d6,
        0x938e397f71ac7853,
        0x44bf383fa903e6af,
        0x164c8ea330587ff8,
    ],
    [
        0x7e2a9673a4a8c4a3,
        0x541007309d5b9410,
        0x5b04ad6376d9cc4a,
        0x06a8d269ce3ed784,
        0x051630c3bff8f84c,
        0x04a3213cdf05ef68,
    ],
    [
        0x9ee764ca2737fe3d,
        0x21567e71e82fc197,
        0x58f94826de4f38af,
        0x45be48e20ffa88fe,
        0x3884331f6c0481da,
        0x0c92415525eec7a8,
    ],
    [
        0x545d1d7a87851206,
        0xebd6351859c5d604,
        0x7a3d74fd7bba4c21,
        0xbbbd3e78ce0bbc17,
        0x1b486acbd1b24bd2,
        0x16381c755aad088d,
    ],
    [
        0x49133e21f2f1a7dc,
        0xf20d881a64afbe3d,
        0x93aa7bc4211c4ec8,
        0x43d1d33d02674d28,
        0x10d3f34f2b511f75,
        0x0773ac208ce0e2b9,
    ],
    [
        0x6bd032334ea9e82b,
        0x8d6af902e96a90e9,
        0xc2d06071c909839b,
        0xbace8929c8ce8fe8,
        0x1308d58c8898e6d4,
        0x0d75cf63b3dda3de,
    ],
    [
        0x72db2a50658d767b,
        0x8abf91faa257b3d5,
        0xe969d6833764ab47,
        0x464170142a1009eb,
        0xb14f01aadb30be2f,
        0x18ae6a856f40715d,
    ],
];

pub const ISO11A_XDEN: [[u64; 6]; 10] = [
    [
        0x334514c97b9c69c2,
        0x40e46332aaa592b3,
        0x1be4e6cb4c145cad,
        0x05007adec6ffdfd2,
        0x7662a93ea0e59c4c,
        0x0cdb4aacf2700470,
    ],
    [
        0xb3aa7a19b3a3f0fc,
        0x1399f1672e854327,
        0x10ae7ada2ec7cac5,
        0x69798c3ff7fdb928,
        0x18e99fa23e6cdb62,
        0x0f242a68d3cf571a,
    ],
    [
        0x39beb66299b9330f,
        0x96a672a1a867c2a0,
        0xd8384c88dbf24b79,
        0x8323143895bc9531,
        0x7f74db3bde4b5777,
        0x011e0724373c9d16,
    ],
    [
        0x7acb03dd30591f6b,
        0x34d5663c660041f8,
        0x257766b848c0fb02,
        0x4b9e0e8bdc65b4be,
        0x465af31c0dd380b2,
        0x026125d7875217f3,
    ],
    [
        0x10b79bc321f54213,
        0xb7b1aee624a0ba95,
        0x0a4af51a1c206aab,
        0xf8ae25babca346df,
        0x0062bd0273e9b169,
        0x13c3f04940437572,
    ],
    [
        0x512a62b3b2f46564,
        0xf0a5f76ebed7a60b,
        0xc4cfbaa922abd93b,
        0xe55616f93494ff6b,
        0x3d811d239fcd5936,
        0x04b92f2fed85ba83,
    ],
    [
        0xe33cea5df1557a19,
        0x2848436246b61323,
        0x97002ab980c93c65,
        0x9488df03d5df8132,
        0xe6161f82703e746f,
        0x1409ba881f35cbc3,
    ],
    [
        0x41a4671d1c9f45b7,
        0x61d4c3abc4c436a3,
        0x3f14b43445246304,
        0xf390f0f8b29102c4,
        0xce70ef04506eda6d,
        0x02731bffa4f75069,
    ],
    [
        0x87d9e0b6d1417b6d,
        0xbc78dd7f5c2e2e8f,
        0x2721d72a3a4c0e85,
        0x329e1b24e47cd4a5,
        0x46f5a8948deecb9a,
        0x0210a9ec42f433a2,
    ],
    [
        0x4b0e194908946a3f,
        0xc8cfde01f1c48bcf,
        0xda7620537f94103f,
        0x447fbbb4945a243b,
        0xb63e62150059bafc,
        0x0df34d04acb09d20,
    ],
];

pub const ISO11A_YNUM: [[u64; 6]; 16] = [
    [
        0x5f0a9f2e39764942,
        0xafd294f9c17c5c16,
        0x0e2a8e2f37ece413,
        0x4b316520df4e3f90,
        0x6ccc0ff2bdffc8ff,
        0x05237cfd8f2b7f90,
    ],
    [
        0xa201b4817b261aa9,
        0x75016e97ad87c543,
        0x67004806eb7d2178,
        0x33a43e9d839ad8a7,
        0xa952d6d1768bdc30,
        0x0ebec265c839f3c3,
    ],
    [
        0x95e3349de1e8f636,
        0x1dd7f4121c2706f3,
        0xceb121d16248cdb2,
        0xcf8de21d8ad874ce,
        0x14d27f914e7ddcbc,
        0x14650f45f353e93e,
    ],
    [
        0xe4d7a55093b4f583,
        0x0a303ff640c800e7,
        0x7570ac4b3d08b4a7,
        0xe3f30d17f6f27efc,
        0x8169456cadb364dc,
        0x0f515fd01bff88d5,
    ],
    [
        0x5e81c6a75b421462,
        0x6297237694089c51,
        0x21e4794fa08ea6ed,
        0xa988c578c95334a3,
        0x1f39f37a6b97dec4,
        0x15d26720b82c67a3,
    ],
    [
        0x7d4481748d755a9c,
        0x2e94e59bf1ded9ed,
        0xda7cd5f4829cfab5,
        0x4d1cfcaf8f997166,
        0x5ee1aa3729abe602,
        0x1591969a2faf474a,
    ],
    [
        0xcba6a9294b2a3c18,
        0x6dd1da0b37e10647,
        0xefcb554541b0a9b6,
        0xf907d0e19a9b548d,
        0x3a1d12447e9b4578,
        0x03f5908784192f5d,
    ],
    [
        0xca38b8e2f076cde4,
        0xe62b682b616bc831,
        0xa71c5b4fbb55538a,
        0x9d16919a06e17a70,
        0xe9359fac06ef137f,
        0x0791b06aa9c0f2ca,
    ],
    [
        0xc630cd91e6e7f51a,
        0xd2ccfa22108cf0f1,
        0x54f9bf13cb94ece5,
        0x847a459d92881cf4,
        0xc95ef5637b23105f,
        0x0271eee0b53625ad,
    ],
    [
        0x9224c4628bba5d30,
        0x727cea24f6d058bc,
        0x7cedaedab4b35d8b,
        0x2ec9930b3d3b617a,
        0xe2846734f0aa3db6,
        0x0a9a5a4846adcbd0,
    ],
    [
        0x04c3256e33fd612f,
        0x35bd86ec7c8a2419,
        0x441ceb8aaa4e542b,
        0x364aea0d29d5812d,
        0xe0b9120cebe2738c,
        0x1916dfa59c9f14a0,
    ],
    [
        0xe5930e1d11fc9a25,
        0x4103744ceebb42f1,
        0xebbeafe54d910a51,
        0x65e8ed247896f893,
        0xc96078acf14ab58a,
        0x07e42ca25b53676a,
    ],
    [
        0x2603f849339b5354,
        0x8a35d40c49a34866,
        0x6fbe890b1ac9cbe9,
        0xac4947b2661267ef,
        0x4464c2291f028513,
        0x118cbef02f1d1b8e,
    ],
    [
        0xc06e1bd471f05f42,
        0x55d76934f6e4a266,
        0x3ecb2a355a6b971e,
        0x92fc4408eee88de3,
        0xf3d8c7d0acfce145,
        0x14c605fb3c0ab306,
    ],
    [
        0x527008ac07e5b6ba,
        0xa28a91e3a9b1d22e,
        0x1b768db3ad13406c,
        0x049c97ab85db0871,
        0x1f637e690647b16f,
        0x1011db4b9d2c8b42,
    ],
    [
        0x872aa6c17d985097,
        0xeecc53161264562a,
        0x07afe37afff55002,
        0x54759078e5be6838,
        0xc4b92d15db8acca8,
        0x106d87d1b51d13b9,
    ],
];

pub const ISO11A_YDEN: [[u64; 6]; 15] = [
    [
        0x5f0a9f2e39764942,
        0xafd294f9c17c5c16,
        0x0e2a8e2f37ece413,
        0x4b316520df4e3f90,
        0x6ccc0ff2bdffc8ff,
        0x05237cfd8f2b7f90,
    ],
    [
        0xae3de771c935c858,
        0x947cda9ec69203c0,
        0xe027e9e9ebe4b515,
        0x40f8a78022d6e1c7,
        0x96b0c36ef60256f6,
        0x181ea1cf92bd65a9,
    ],
    [
        0x5a6277f1b8e42fe6,
        0xd457d18195d1e061,
        0xa2619aa1b23470a4,
        0x361fbbffa948ff0f,
        0xefc3a53525acea55,
        0x14c246d4e707b036,
    ],
    [
        0x98b8901e03405e42,
        0x085e8173f1f43a7e,
        0x2a339c5fbbc948e8,
        0x0ed2aede4e34b0a8,
        0xc8712f97101ca1ed,
        0x08860b2a3172a8a5,
    ],
    [
        0xe2bd8c6b0eea38fe,
        0xc1e3ab7d2e9970a4,
        0xfa079f581ebeffa8,
        0x064f91218acd1334,
        0xd002bd5d7d66e909,
        0x07a59ebac6331df0,
    ],
    [
        0x7eb6e38f9c6587af,
        0x743f49ae7ecd1a63,
        0xa8c9aed956b448dc,
        0xa166b4a25bd0a57b,
        0x55a1a6702e2642c4,
        0x193a5588b7c84086,
    ],
    [
        0xf172ba83e1d14527,
        0x2416af07182d9be8,
        0xde800a1509962d73,
        0x35f63242e8e2e1df,
        0x64ba56c8725aad7c,
        0x1736897c7f5fb1ba,
    ],
    [
        0x40e6e0c9691297bf,
        0x2239cda512476898,
        0x153ee2c6a38813dd,
        0x499a4071e15cd5a9,
        0xf2600bffd9ceb397,
        0x03840846fbb6a2d1,
    ],
    [
        0xdb80ed873870fc68,
        0xb21e1c3a5ddcc352,
        0x32cdf95c03b0c6a5,
        0xf3447b192f328e82,
        0x8f7d1342f176ebdc,
        0x11577e1b71680e3b,
    ],
    [
        0x441c2436c88daea5,
        0xcfa972cdc91f8d96,
        0xee7a730027de7ae6,
        0x9dcbfab18fbcc693,
        0xa8d3fdd71277184a,
        0x09a580fdf3b07165,
    ],
    [
        0x0761de0643ab871a,
        0x1e3df939b9ceb7a1,
        0x9dc03f53e5b78dd6,
        0xe85d7efd9b10d473,
        0xfc5c0555620e906a,
        0x1632fe8b5eb46fc7,
    ],
    [
        0x70db617692459c70,
        0xafa255ca839a45ab,
        0x726368d2219eb91f,
        0x2708762ba637fbae,
        0x1c16880e43c62f9f,
        0x145a10d09532f4b5,
    ],
    [
        0x464963120a0df36f,
        0xbbccf5287786723f,
        0x4ae40346bdbe2e6a,
        0xacdc24ff39bec7da,
        0xe3e1bd83e33e0499,
        0x03dfde52fce8c5ba,
    ],
    [
        0xc8d2161d8eaafe14,
        0xcc2eca6d37261670,
        0x52e224fffe5ca08f,
        0x50ee0bb72bab3455,
        0x05eed719546aed65,
        0x062600b423d7e9c4,
    ],
    [
        0xdadeab10f8f88bd5,
        0x3faa4c9daeaf18d6,
        0xc5e5b8bcffb1bdbd,
        0xbf1c00f59b33e5c1,
        0x1c726dedbceefa17,
        0x094a5f563d6a6ff4,
    ],
];
