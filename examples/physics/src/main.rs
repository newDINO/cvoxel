use bevy::{prelude::*, render::{mesh::{Indices, PrimitiveTopology, VertexAttributeValues}, render_asset::RenderAssetUsages}};
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use cvoxel::CVoxels;
use nalgebra::{Isometry3, Vector3};
use pvoxel::{PVoxels, RigidType};

fn main() {
    App::new()
        .add_plugins(PanOrbitCameraPlugin)
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .add_systems(Update, update_from_pvoxel_transform)
        .add_systems(Update, update_physics)
        .run();
}

// systems
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // camera
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(0.0, 10.0, 15.0)
                .looking_at(Vec3::new(0.0, 0.0, 0.0), Dir3::Y),
            ..Default::default()
        },
        PanOrbitCamera::default(),
    ));

    // light
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 1e7,
            color: Color::WHITE,
            shadows_enabled: true,
            ..Default::default()
        },
        transform: Transform::from_xyz(0.0, 10.0, 5.0),
        ..Default::default()
    });

    // meshes
    let shapes = [
        Cuboid::new(10.0, 1.0,10.0).mesh().build(),
        Torus::new(1.0, 3.0).mesh().build(),
    ];
    let rigid_tys = [
        RigidType::Fixed,
        RigidType::Dynamic,
    ];
    let poses = [
        Vector3::new(0.0, -0.5, 0.0),
        Vector3::new(0.0, 2.0, 0.0),
    ];
    let material_handle = materials.add(Color::WHITE);
    let mut objs = Vec::with_capacity(shapes.len());
    let mut ids = Vec::with_capacity(shapes.len());
    for i in 0..shapes.len() {
        let mesh = &shapes[i];
        let mut cvoxels = voxelize_mesh(&mesh, 0.1).unwrap();
        cvoxels.transform.translation.vector = poses[i];
        let surface_mesh = cvoxel_surface_mesh(&cvoxels);
        let mesh_handle = meshes.add(surface_mesh);
        let id = commands.spawn(PbrBundle {
            mesh: mesh_handle,
            material: material_handle.clone(),
            ..Default::default()
        }).id();
        objs.push(PVoxels::from_cvoxels(cvoxels, 1.0, rigid_tys[i]));
        ids.push(id);
    }
    commands.insert_resource(Objects {
        objs, ids
    });
}

fn update_from_pvoxel_transform(objects: Res<Objects>, mut query: Query<&mut Transform>) {
    for i in 0..objects.objs.len() {
        *query.get_mut(objects.ids[i]).unwrap() = isometry_to_transform(&objects.objs[i].transform);
    }
}

fn update_physics(objects: ResMut<Objects>, time: Res<Time>) {
    let objects = objects.into_inner();
    let dt = time.delta_seconds();
    // gravity
    for objects in &mut objects.objs {
        if objects.ty == RigidType::Fixed {
            continue;
        }
        objects.vel.y -= dt * 9.8;
    }

    // collision handling

    // time integral
    for objects in &mut objects.objs {
        if objects.ty == RigidType::Fixed {
            continue;
        }
        objects.step_dt(dt);
    }
}

// structs
#[derive(Resource)]
struct Objects {
    objs: Vec<PVoxels>,
    ids: Vec<Entity>
}

// utils
fn isometry_to_transform(
    isometry: &Isometry3<f32>,
) -> Transform {
    let mut transform = Transform::IDENTITY;
    transform.translation = Vec3::from_array(isometry.translation.vector.data.0[0]);
    transform.rotation = Quat::from_array(isometry.rotation.quaternion().coords.data.0[0]);
    transform
}

fn voxelize_mesh(mesh: &Mesh, dx: f32) -> Option<CVoxels> {
    let mesh_attr = mesh.attribute(Mesh::ATTRIBUTE_POSITION).unwrap();
    if let VertexAttributeValues::Float32x3(v) = mesh_attr {
        if let Some(indices) = mesh.indices() {
            match indices {
                Indices::U16(ids) => CVoxels::from_indexed_mesh(&v, &ids, dx),
                Indices::U32(ids) => CVoxels::from_indexed_mesh(&v, &ids, dx),
            }
        } else {
            CVoxels::from_trimesh(&v, dx)
        }
    } else {
        None
    }
}

fn cvoxel_surface_mesh(voxels: &CVoxels) -> Mesh {
    let surface_mesh = voxels.surface_mesh();
    Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    )
    .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, surface_mesh.position)
    .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, surface_mesh.color)
    .with_computed_normals()
}