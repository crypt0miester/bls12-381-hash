// BLS12-381 G2 SSWU + iso-3 + psi constants, Montgomery form, limbs
// little-endian; Fp2 as [c0, c1]. Extracted from bls12_381 0.8.0.

pub const SSWU2_ELLP_A: [[u64; 6]; 2] = [
    [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
    [0xe53a000003135242, 0x01080c0fdef80285, 0xe7889edbe340f6bd, 0x0b51375126310601, 0x02d6985717c744ab, 0x1220b4e979ea5467],
];

pub const SSWU2_ELLP_B: [[u64; 6]; 2] = [
    [0x22ea00000cf89db2, 0x6ec832df71380aa4, 0x6e1b94403db5a66e, 0x75bf3c53a79473ba, 0x3dd3a569412c0a34, 0x125cdb5e74dc4fd1],
    [0x22ea00000cf89db2, 0x6ec832df71380aa4, 0x6e1b94403db5a66e, 0x75bf3c53a79473ba, 0x3dd3a569412c0a34, 0x125cdb5e74dc4fd1],
];

pub const SSWU2_XI: [[u64; 6]; 2] = [
    [0x87ebfffffff9555c, 0x656fffe5da8ffffa, 0x0fd0749345d33ad2, 0xd951e663066576f4, 0xde291a3d41e980d3, 0x0815664c7dfe040d],
    [0x43f5fffffffcaaae, 0x32b7fff2ed47fffd, 0x07e83a49a2e99d69, 0xeca8f3318332bb7a, 0xef148d1ea0f4c069, 0x040ab3263eff0206],
];

pub const SSWU2_C1_NEG_B_OVER_A: [[u64; 6]; 2] = [
    [0x903c555555474fb3, 0x5f98cc95ce451105, 0x9f8e582eefe0fade, 0xc68946b6aebbd062, 0x467a4ad10ee6de53, 0x0e7146f483e23a05],
    [0x29c2aaaaaab85af8, 0xbf133368e30eeefa, 0xc7a27a7206cffb45, 0x9dee04ce44c9425c, 0x04a15ce53464ce83, 0x0b8fcaf5b59dac95],
];

pub const ISO3_XNUM: [[[u64; 6]; 2]; 4] = [
    [
        [0x47f671c71ce05e62, 0x06dd57071206393e, 0x7c80cd2af3fd71a2, 0x048103ea9e6cd062, 0xc54516acc8d037f6, 0x13808f550920ea41],
        [0x47f671c71ce05e62, 0x06dd57071206393e, 0x7c80cd2af3fd71a2, 0x048103ea9e6cd062, 0xc54516acc8d037f6, 0x13808f550920ea41],
    ],
    [
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x5fe55555554c71d0, 0x873fffdd236aaaa3, 0x6a6b4619b26ef918, 0x21c2888408874945, 0x2836cda7028cabc5, 0x0ac73310a7fd5abd],
    ],
    [
        [0x0a0c5555555971c3, 0xdb0c00101f9eaaae, 0xb1fb2f941d797997, 0xd3960742ef416e1c, 0xb70040e2c20556f4, 0x149d7861e581393b],
        [0xaff2aaaaaaa638e8, 0x439fffee91b55551, 0xb535a30cd9377c8c, 0x90e144420443a4a2, 0x941b66d3814655e2, 0x0563998853fead5e],
    ],
    [
        [0x40aac71c71c725ed, 0x190955557a84e38e, 0xd817050a8f41abc3, 0xd86485d4c87f6fb1, 0x696eb479f885d059, 0x198e1a74328002d2],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
    ],
];

pub const ISO3_XDEN: [[[u64; 6]; 2]; 3] = [
    [
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x1f3affffff13ab97, 0xf25bfc611da3ff3e, 0xca3757cb3819b208, 0x3e6427366f8cec18, 0x03977bc86095b089, 0x04f69db13f39a952],
    ],
    [
        [0x447600000027552e, 0xdcb8009a43480020, 0x6f7ee9ce4a6e8b59, 0xb10330b7c0a95bc6, 0x6140b1fcfb1e54b7, 0x0381be097f0bb4e1],
        [0x7588ffffffd8557d, 0x41f3ff646e0bffdf, 0xf7b1e8d2ac426aca, 0xb3741acd32dbb6f8, 0xe9daf5b9482d581f, 0x167f53e0ba7431b8],
    ],
    [
        [0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba, 0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
    ],
];

pub const ISO3_YNUM: [[[u64; 6]; 2]; 4] = [
    [
        [0x96d8f684bdfc77be, 0xb530e4f43b66d0e2, 0x184a88ff379652fd, 0x57cb23ecfae804e1, 0x0fd2e39eada3eba9, 0x08c8055e31c5d5c3],
        [0x96d8f684bdfc77be, 0xb530e4f43b66d0e2, 0x184a88ff379652fd, 0x57cb23ecfae804e1, 0x0fd2e39eada3eba9, 0x08c8055e31c5d5c3],
    ],
    [
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0xbf0a71c71c91b406, 0x4d6d55d28b7638fd, 0x9d82f98e5f205aee, 0xa27aa27b1d1a18d5, 0x02c3b2b2d2938e86, 0x0c7d13420b09807f],
    ],
    [
        [0xd7f9555555531c74, 0x21cffff748daaaa8, 0x5a9ad1866c9bbe46, 0x4870a2210221d251, 0x4a0db369c0a32af1, 0x02b1ccc429ff56af],
        [0xe205aaaaaaac8e37, 0xfcdc000768795556, 0x0c96011a8a1537dd, 0x1c06a963f163406e, 0x010df44c82a881e6, 0x174f45260f808feb],
    ],
    [
        [0xa470bda12f67f35c, 0xc0fe38e23327b425, 0xc9d3d0f2c6f0678d, 0x1c55c9935b5a982e, 0x27f6c0e2f0746764, 0x117c5e6e28aa9054],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
    ],
];

pub const ISO3_YDEN: [[[u64; 6]; 2]; 4] = [
    [
        [0x0162fffffa765adf, 0x8f7bea480083fb75, 0x561b3c2259e93611, 0x11e19fc1a9c875d5, 0xca713efc00367660, 0x03c6a03d41da1151],
        [0x0162fffffa765adf, 0x8f7bea480083fb75, 0x561b3c2259e93611, 0x11e19fc1a9c875d5, 0xca713efc00367660, 0x03c6a03d41da1151],
    ],
    [
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x5db0fffffd3b02c5, 0xd713f52358ebfdba, 0x5ea60761a84d161a, 0xbb2c75a34ea6c44a, 0x0ac6735921c1119b, 0x0ee3d913bdacfbf6],
    ],
    [
        [0x66b10000003affc5, 0xcb1400e764ec0030, 0xa73e5eb56fa5d106, 0x8984c913a0fe09a9, 0x11e10afb78ad7f13, 0x05429d0e3e918f52],
        [0x534dffffffc4aae6, 0x5397ff174c67ffcf, 0xbff273eb870b251d, 0xdaf2827152870915, 0x393a9cbaca9e2dc3, 0x14be74dbfaee5748],
    ],
    [
        [0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba, 0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
    ],
];

pub const PSI_X_C1: [u64; 6] = [
    0x890dc9e4867545c3,
    0x2af322533285a5d5,
    0x50880866309b7e2c,
    0xa20d1b8c7e881024,
    0x14e4f04fe2db9068,
    0x14e56d3f1564853a,
];

pub const PSI_Y: [[u64; 6]; 2] = [
    [0x3e2f585da55c9ad1, 0x4294213d86c18183, 0x382844c88b623732, 0x92ad2afd19103e18, 0x1d794e4fac7cf0b9, 0x0bd592fc7d825ec8],
    [0x7bcfa7a25aa30fda, 0xdc17dec12a927e7c, 0x2f088dd86b4ebef1, 0xd1ca2087da74d4a7, 0x2da2596696cebc1d, 0x0e2b7eedbbfd87d2],
];

pub const PSI2_X_C0: [u64; 6] = [
    0xcd03c9e48671f071,
    0x5dab22461fcda5d2,
    0x587042afd3851b95,
    0x8eb60ebe01bacb9e,
    0x03f97d6e83d050d2,
    0x18f0206554638741,
];

pub const C256_MONT: [u64; 6] = [
    0x075b3cd7c5ce820f,
    0x3ec6ba621c3edb0b,
    0x168a13d82bff6bce,
    0x87663c4bf8c449d2,
    0x15f34c83ddc8d830,
    0x0f9628b49caa2e85,
];


// Shallue-van de Woestijne direct-map constants for E2 (Wahby-Boneh
// eprint 2019/403 section 3, u0 = -1), Montgomery form, derived offline.
pub const SVDW2_F_U0: [[u64; 6]; 2] = [
    [
        0xee1d00000009aaa1,
        0x86840025e97c0007,
        0x4f7823c40df41de8,
        0x9e7c71f069ece051,
        0x7dde005a606d6b99,
        0x0de0f8777c82e085,
    ],
    [
        0xaa270000000cfff3,
        0x53cc0032fc34000a,
        0x478fe97a6b0a807f,
        0xb1d37ebee6ba24d7,
        0x8ec9733bbf78ab2f,
        0x09d645513d83de7e,
    ],
];

pub const SVDW2_B: [[u64; 6]; 2] = [
    [
        0xaa270000000cfff3,
        0x53cc0032fc34000a,
        0x478fe97a6b0a807f,
        0xb1d37ebee6ba24d7,
        0x8ec9733bbf78ab2f,
        0x09d645513d83de7e,
    ],
    [
        0xaa270000000cfff3,
        0x53cc0032fc34000a,
        0x478fe97a6b0a807f,
        0xb1d37ebee6ba24d7,
        0x8ec9733bbf78ab2f,
        0x09d645513d83de7e,
    ],
];

pub const SVDW2_SQRT_M3: [u64; 6] = [
    0x1dec6c36f3181f22,
    0xb4b9bb641054b457,
    0x25695a2be9415286,
    0x982b6cbf66c749bc,
    0x7d58e1ae1feb7873,
    0x062c96300937c0b9,
];

pub const SVDW2_C1: [u64; 6] = [
    0xecfb361b798dba3a,
    0xc100ddb891865a2c,
    0x0ec08ff1232bda8e,
    0xd5c13cc6f1ca4721,
    0x47222a47bf7b5c04,
    0x0110f184e51c5f59,
];

pub const SVDW2_INV_3U0SQ: [u64; 6] = [
    0x4e02555555561c71,
    0x0dc400030ce6aaab,
    0xb9e369ddc0631701,
    0xc03efa7472742996,
    0xa614ce0162fa175e,
    0x18a82b8824803b42,
];

// Adapted iso-3 evaluation constants: each cubic is the identity
// (y + g)(y^2 + q1) + q0 (plus the leading coefficient when not
// monic); derived and checked offline.

pub const ISO3A_XNUM: [[[u64; 6]; 2]; 4] = [
    [
        [
            0x447600000027552e,
            0xdcb8009a43480020,
            0x6f7ee9ce4a6e8b59,
            0xb10330b7c0a95bc6,
            0x6140b1fcfb1e54b7,
            0x0381be097f0bb4e1,
        ],
        [
            0x7588ffffffd8557d,
            0x41f3ff646e0bffdf,
            0xf7b1e8d2ac426aca,
            0xb3741acd32dbb6f8,
            0xe9daf5b9482d581f,
            0x167f53e0ba7431b8,
        ],
    ],
    [
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
        [
            0x3112ffffffb1004f,
            0x653bfeca2ac3ffbf,
            0x8832ff0461d3df70,
            0x0270ea1572325b32,
            0x889a43bc4d0f0368,
            0x12fd95d73b687cd7,
        ],
    ],
    [
        [
            0x68c3000007964dbf,
            0xdafc1dc1b5040639,
            0x33b5ba30e20d67d9,
            0x38a40ccd12064556,
            0x698596623c80d54c,
            0x19a539a535c21604,
        ],
        [
            0x68c3000007964dbf,
            0xdafc1dc1b5040639,
            0x33b5ba30e20d67d9,
            0x38a40ccd12064556,
            0x698596623c80d54c,
            0x19a539a535c21604,
        ],
    ],
    [
        [
            0x40aac71c71c725ed,
            0x190955557a84e38e,
            0xd817050a8f41abc3,
            0xd86485d4c87f6fb1,
            0x696eb479f885d059,
            0x198e1a74328002d2,
        ],
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
    ],
];

pub const ISO3A_YNUM: [[[u64; 6]; 2]; 4] = [
    [
        [
            0x66b10000003affc5,
            0xcb1400e764ec0030,
            0xa73e5eb56fa5d106,
            0x8984c913a0fe09a9,
            0x11e10afb78ad7f13,
            0x05429d0e3e918f52,
        ],
        [
            0x534dffffffc4aae6,
            0x5397ff174c67ffcf,
            0xbff273eb870b251d,
            0xdaf2827152870915,
            0x393a9cbaca9e2dc3,
            0x14be74dbfaee5748,
        ],
    ],
    [
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
        [
            0x4bd8fffffc9dae0d,
            0x6433f2ba4bcbfd39,
            0xa0aa60287e92e8b3,
            0xf71fb2c44c015530,
            0x85c3ab653547bebc,
            0x00dce0edc17e2870,
        ],
    ],
    [
        [
            0x05d200003345ccba,
            0x0ae8c91755182a10,
            0x2c110485dfba800a,
            0x4a019765b4a23a51,
            0xdc147e44919f7c44,
            0x076a517f64429e8a,
        ],
        [
            0x05d200003345ccba,
            0x0ae8c91755182a10,
            0x2c110485dfba800a,
            0x4a019765b4a23a51,
            0xdc147e44919f7c44,
            0x076a517f64429e8a,
        ],
    ],
    [
        [
            0xa470bda12f67f35c,
            0xc0fe38e23327b425,
            0xc9d3d0f2c6f0678d,
            0x1c55c9935b5a982e,
            0x27f6c0e2f0746764,
            0x117c5e6e28aa9054,
        ],
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
    ],
];

pub const ISO3A_YDEN: [[[u64; 6]; 2]; 3] = [
    [
        [
            0x66b10000003affc5,
            0xcb1400e764ec0030,
            0xa73e5eb56fa5d106,
            0x8984c913a0fe09a9,
            0x11e10afb78ad7f13,
            0x05429d0e3e918f52,
        ],
        [
            0x534dffffffc4aae6,
            0x5397ff174c67ffcf,
            0xbff273eb870b251d,
            0xdaf2827152870915,
            0x393a9cbaca9e2dc3,
            0x14be74dbfaee5748,
        ],
    ],
    [
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
        [
            0x5db0fffffd3b02c5,
            0xd713f52358ebfdba,
            0x5ea60761a84d161a,
            0xbb2c75a34ea6c44a,
            0x0ac6735921c1119b,
            0x0ee3d913bdacfbf6,
        ],
    ],
    [
        [
            0x68e600002c4c7e5e,
            0xc178adbd5e882457,
            0x1d87c42f1e183bbb,
            0x39e198fc98c676d4,
            0x42ad578c84e3a6ae,
            0x15cd21ea642f42a6,
        ],
        [
            0x68e600002c4c7e5e,
            0xc178adbd5e882457,
            0x1d87c42f1e183bbb,
            0x39e198fc98c676d4,
            0x42ad578c84e3a6ae,
            0x15cd21ea642f42a6,
        ],
    ],
];

pub const ISO3A_XDEN: [[[u64; 6]; 2]; 2] = [
    [
        [
            0x447600000027552e,
            0xdcb8009a43480020,
            0x6f7ee9ce4a6e8b59,
            0xb10330b7c0a95bc6,
            0x6140b1fcfb1e54b7,
            0x0381be097f0bb4e1,
        ],
        [
            0x7588ffffffd8557d,
            0x41f3ff646e0bffdf,
            0xf7b1e8d2ac426aca,
            0xb3741acd32dbb6f8,
            0xe9daf5b9482d581f,
            0x167f53e0ba7431b8,
        ],
    ],
    [
        [
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ],
        [
            0x1f3affffff13ab97,
            0xf25bfc611da3ff3e,
            0xca3757cb3819b208,
            0x3e6427366f8cec18,
            0x03977bc86095b089,
            0x04f69db13f39a952,
        ],
    ],
];
